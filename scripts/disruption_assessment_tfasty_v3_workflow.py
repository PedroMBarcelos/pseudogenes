#!/usr/bin/env python3
"""Version 3: cross-checked TFASTY plus unmasked genomic terminal-codon assessment."""
import argparse, json, re, subprocess, tempfile
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import os
from Bio.Seq import Seq
from Bio.Data import CodonTable

FIELDS = ('frameshifts','premature_stop_codons','missing_start_codon','missing_stop_codon')

def parse_m10(text):
    hits=[]; hit=None; block=None
    for line in text.splitlines():
        if line.startswith('>>>'):
            hit=block=None
        elif line.startswith('>>'):
            hit={'id':line[2:].strip().split()[0], 'meta':{}, 'blocks':[]}; hits.append(hit); block=None
        elif line.strip() == '>--':
            if hit is None: raise ValueError('secondary alignment without target')
            hit={'id':hit['id'], 'meta':{}, 'blocks':[]}; hits.append(hit); block=None
        elif hit is not None and line.startswith('>'):
            block={'meta':{},'seq':''}; hit['blocks'].append(block)
        elif hit is not None and line.startswith(';'):
            key, sep, value=line[1:].partition(':'); key=key.strip()
            if key=='al_cons': block=None
            elif sep: (block['meta'] if block is not None else hit['meta'])[key]=value.strip()
        elif block is not None and re.fullmatch(r'[A-Za-z*./\\-]+',line.strip()):
            block['seq']+=line.strip().upper()
    return hits

def inspect(hit, protein, dna):
    q,s=hit['blocks']; qm,sm=q['meta'],s['meta']; qa,sa=q['seq'],s['seq']
    qs,qe,ss,se=map(int,(qm['al_start'],qm['al_stop'],sm['al_start'],sm['al_stop']))
    if hit['meta'].get('tfy_frame')!='f' or not (1<=qs<=qe<=len(protein) and 1<=ss<=se<=len(dna)):
        raise ValueError('unexpected alignment orientation or bounds')
    if int(qm['al_display_start'])!=qs or int(sm['al_display_start'])!=ss or len(qa)!=len(sa) or qa.replace('-','')!=protein[qs-1:qe]:
        raise ValueError('alignment display includes flanks or inconsistent query')
    if not re.fullmatch(r'[A-Z*\\/\-]+',sa) or not all(a=='-' for a,b in zip(qa,sa) if b in '/\\'):
        raise ValueError('unsupported traceback symbols')
    full=qs==1 and qe==len(protein)
    fs=sa.count('/')+sa.count('\\'); stops=sum(b=='*' and a not in '-*' for a,b in zip(qa,sa))
    # Positive calls are traceback evidence, not a definitive biological diagnosis.
    vals=dict(frameshifts=fs if fs or full else None,premature_stop_codons=stops if stops or (full and '*' not in sa) else None,missing_start_codon=None,missing_stop_codon=None)
    # Terminal calls only on a verified ungapped translation. The benchmark
    # demonstrates al_stop is three bases before the last codon start (1-based).
    if not fs and '-' not in qa+sa:
        a=ss-1; b=a+3*len(sa); coding=dna[a:b]
        if b<=len(dna) and set(coding)<=set('ACGT') and str(Seq(coding).translate(table=11))==sa and se==b-3:
            ct=CodonTable.unambiguous_dna_by_id[11]
            if qs==1 and coding[:3] in ct.start_codons: vals['missing_start_codon']=0
            if qe==len(protein) and dna[b:b+3] in ct.stop_codons: vals['missing_stop_codon']=0
    return vals, {'query_start':qs,'query_end':qe,'target_start':ss,'target_stop_reported':se,'full_reference':full,'query_alignment':qa,'target_alignment':sa,'alignment_metadata':hit['meta'],'scope':'local traceback; positive calls are candidate evidence; terminal absence not called'}

def load_refs(path, wanted):
    refs={}; name=None; seq=[]
    def flush():
        if name is not None and name in wanted:
            if name in refs: raise ValueError('duplicate FASTA identifier '+name)
            refs[name]=''.join(seq).upper()
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                flush(); name=line[1:].split()[0]; seq=[]
            elif name in wanted: seq.append(line.strip())
        flush()
    return refs

def parse_m8(text):
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) != 13:
            raise ValueError('unexpected m8CC column count: '+str(len(f)))
        rows.append((f[1], tuple(map(int, f[6:10])), f[12]))
    return rows


def traceback_code(qa, sa):
    ops = []
    for q, s in zip(qa, sa):
        if s in '/\\': op = 'F' if s == '/' else 'R'
        elif s == '-' and q != '-': op = 'D'
        elif q == '-' and s != '-': op = 'I'
        elif q != '-' and s != '-': op = 'M'
        else: raise ValueError('double gap')
        if ops and ops[-1][1] == op: ops[-1][0] += 1
        else: ops.append([1, op])
    return ''.join(f'{n}{op}' for n, op in ops)


def verify(hit, m8, protein, dna, region):
    vals, detail = inspect(hit, protein, dna)
    qa, sa = detail['query_alignment'], detail['target_alignment']
    code = traceback_code(qa, sa)
    coords = tuple(detail[k] for k in ('query_start','query_end','target_start','target_stop_reported'))
    if m8 != (hit['id'], coords, code):
        raise ValueError('m8CC and m10 disagree: '+str(m8)+' vs '+str((hit['id'],coords,code)))
    start, end = int(region['start']), int(region['end'])
    strand = region['strand']
    if strand not in ('+','-') or start < 0 or end-start != len(dna):
        raise ValueError('invalid genomic interval or strand')
    def genomic(a,b):
        return [start+a,start+b] if strand=='+' else [end-b,end-a]
    a = detail['target_start']-1
    qpos = detail['query_start']-1
    shifted = False; events = []; stops = []; unchecked = []; uncertain_stop = False; ambiguous_dna = False
    first = last = None
    for i,(q,s) in enumerate(zip(qa,sa)):
        if q != '-': qpos += 1
        if s in '/\\':
            old=a; a += -1 if s=='/' else 1
            if not 0<=a<=len(dna): raise ValueError('frame transition outside candidate')
            events.append({'operation':'F' if s=='/' else 'R', 'after_reference_aa':qpos,
                           'oriented_cursor_before':old,'oriented_cursor_after':a,
                           'genomic_transition_interval_bed':genomic(min(old,a),max(old,a)),
                           'localization':'traceback cursor; not exact indel coordinate'})
            shifted=True; continue
        if s=='-': continue
        codon=dna[a:a+3]
        if a<0 or len(codon)!=3: raise ValueError('codon outside candidate')
        ambiguous_dna |= not set(codon)<=set('ACGT')
        translated=str(Seq(codon).translate(table=11))
        if shifted:
            unchecked.append(i+1)
        elif translated!=s:
            raise ValueError(f'DNA translation mismatch at alignment column {i+1}: {codon} -> {translated}, expected {s}')
        if q!='-':
            if first is None: first=(qpos,a,codon,shifted)
            last=(qpos,a,codon,shifted)
        if s=='*':
            if q not in '-*' and not shifted and translated=='*':
                stops.append({'reference_aa':qpos,'codon':codon,'oriented_start0':a,
                              'genomic_start0':genomic(a,a+3)[0], 'genomic_end0':genomic(a,a+3)[1],
                              'chromosome':region.get('chromosome'),'strand':strand})
            else: uncertain_stop=True
        a+=3; shifted=False
    if shifted or a-3!=detail['target_stop_reported']:
        raise ValueError('translated target endpoint inconsistent with traceback')
    unique={(x['genomic_start0'],x['genomic_end0']):x for x in stops}
    full=detail['full_reference']; fs=len(events)
    vals['frameshifts']=fs if fs or (full and not ambiguous_dna) else None
    vals['premature_stop_codons']=len(unique) if unique else (0 if full and not uncertain_stop and not ambiguous_dna else None)
    # Positive terminal absence remains uncalled; verify presence at aligned endpoints only.
    vals['missing_start_codon']=None; vals['missing_stop_codon']=None
    ct=CodonTable.unambiguous_dna_by_id[11]
    if first and qa[0]!='-' and sa[0] not in '-/\\' and first[0]==1 and not first[3] and first[2] in ct.start_codons:
        vals['missing_start_codon']=0
    if last and qa[-1]!='-' and sa[-1] not in '-/\\' and last[0]==len(protein) and not last[3] and dna[a:a+3] in ct.stop_codons:
        vals['missing_stop_codon']=0
    detail.update({'alignment_code':code,'reference_coverage_pct':100*(coords[1]-coords[0]+1)/len(protein),
                   'frameshift_events':events,'verified_stops':list(unique.values()),'unverified_stop_symbol':uncertain_stop,'ambiguous_aligned_dna':ambiguous_dna,
                   'unchecked_frame_transition_columns':unchecked,'target_end0':a,
                   'scope':'local traceback; frame transitions cross-checked, exact indel positions unresolved; terminal absence not called'})
    return vals,detail


def select_primary(hits, rows, protein, dna):
    """Match first reported m10 traceback, never aggregate local alignments."""
    if not hits or hits[0]['id'] != 'target':
        raise ValueError('no primary target alignment')
    hit = hits[0]
    _, detail = inspect(hit, protein, dna)
    coords = tuple(detail[k] for k in ('query_start','query_end','target_start','target_stop_reported'))
    code = traceback_code(detail['query_alignment'], detail['target_alignment'])
    key = (hit['id'], coords, code)
    matches = [row for row in rows if row == key]
    if len(matches) != 1:
        raise ValueError(f'primary traceback has {len(matches)} matching m8CC rows')
    return hit, matches[0]


def assess_boundaries(r, genomes, window):
    """BED coordinates; offsets relative to strand-oriented candidate DNA.

    Noncanonical expected sites are evidence, not proof of terminal absence.
    Alternative codons are descriptive and never used to assert functionality.
    """
    e = r.setdefault('mutation_evidence', {})
    t = e.get('tfasty', {})
    b = {'status': 'not_assessed', 'genetic_code': 11,
         'alternative_window_nt': window, 'start': {}, 'stop': {}}
    e['terminal_boundary'] = b
    if 'query_alignment' not in t:
        b['reason'] = 'alignment_unavailable'; return
    try:
        region = r['region']; chrom = region['chromosome']
        lo, hi = int(region['start']), int(region['end'])
        strand = region['strand']; dna = r['sequence_nt'].upper()
        seq = genomes.get(chrom)
        if seq is None: raise ValueError('contig_not_found')
        if not (0 <= lo < hi <= len(seq)) or strand not in ('+', '-'):
            raise ValueError('invalid_genomic_interval')
        def fetch(offset):
            x, y = (lo+offset, lo+offset+3) if strand == '+' else (hi-offset-3, hi-offset)
            if x < 0 or y > len(seq): return None
            codon = seq[x:y]
            if strand == '-': codon = str(Seq(codon).reverse_complement())
            return dict(codon=codon, genomic_start0=x, genomic_end0=y,
                        chromosome=chrom, strand=strand, oriented_start0=offset,
                        outside_candidate=(offset < 0 or offset+3 > len(dna)))
        actual = seq[lo:hi] if strand == '+' else str(Seq(seq[lo:hi]).reverse_complement())
        if actual != dna: raise ValueError('candidate_genome_sequence_mismatch')
        b['status'] = 'genome_verified'
        qa, sa = t['query_alignment'], t['target_alignment']
        ct = CodonTable.unambiguous_dna_by_id[11]
        unchecked = set(t.get('unchecked_frame_transition_columns', []))
        specs = [
            ('start', t['target_start']-1, t['query_start']==1 and qa[0]!='-' and sa[0] not in '-/\\' and 1 not in unchecked, ct.start_codons),
            ('stop', t['target_end0'], t['query_end']==len_reference(r) and qa[-1]!='-' and sa[-1] not in '-/\\' and len(sa) not in unchecked, ct.stop_codons)
        ]
        for kind, offset, anchored, valid in specs:
            d = b[kind] = {'expected_site_status':'unknown', 'expected_site_altered':None, 'alternatives':[]}
            if not anchored:
                d['reason']='reference_terminus_not_reliably_anchored'; continue
            site = fetch(offset)
            if site is None:
                d['reason']='contig_boundary'; continue
            d['expected_site']=site
            if not set(site['codon']) <= set('ACGT'):
                d['reason']='ambiguous_codon'; continue
            present = site['codon'] in valid
            d['expected_site_status']='canonical' if present else 'noncanonical'
            d['expected_site_altered']=int(not present)
            if present:
                r['missing_'+kind+'_codon']=0
                continue
            # Scan in the local terminal frame only. Stop at a barrier or
            # unavailable/ambiguous DNA; do not cross an intervening stop.
            for direction in (-1, 1):
                for delta in range(3, window+1, 3):
                    alt=fetch(offset+direction*delta)
                    if alt is None or not set(alt['codon']) <= set('ACGT'): break
                    if alt['codon'] in valid:
                        alt['offset_from_expected_nt']=direction*delta
                        d['alternatives'].append(alt)
                        if kind=='stop': break
                    if alt['codon'] in ct.stop_codons: break
            d['reason']='expected_site_altered; biological_terminal_absence_not_inferred'
        b['status']='assessed' if any(b[k].get('expected_site_altered') is not None for k in ('start','stop')) else 'no_anchored_site'
    except (KeyError, ValueError, TypeError, IndexError) as exc:
        b['status']='genome_validation_failed'; b['reason']=str(exc)


def len_reference(r):
    lengths={int(h['qlen']) for h in r['best_protein']['hsps']}
    if len(lengths)!=1: raise ValueError('inconsistent_reference_length')
    return lengths.pop()


def run_record(r, refs, a):
    ev=r.setdefault('mutation_evidence',{})
    stage='input_validation'
    try:
        pid=r['best_protein']['protein_id']; dna=r['sequence_nt'].upper()
        if pid not in refs: raise ValueError('reference not found: '+pid)
        protein=refs[pid]
        if any(int(h['qlen'])!=len(protein) for h in r['best_protein']['hsps']): raise ValueError('reference length mismatch')
        if len(dna)!=int(r['region']['end'])-int(r['region']['start']): raise ValueError('candidate length mismatch')
        with tempfile.TemporaryDirectory() as td:
            q=Path(td)/'query.faa'; s=Path(td)/'target.fna'
            q.write_text('>query\n'+protein+'\n'); s.write_text('>target\n'+dna+'\n')
            outputs={}
            for mode in ('10','8CC'):
                stage='execution_m'+mode
                run=subprocess.run([str(a.tfasty.resolve()),'-3','-m',mode,str(q),str(s)],capture_output=True,text=True,check=True,timeout=a.timeout)
                outputs[mode]=run.stdout
            stage='parse_m10'; hits=parse_m10(outputs['10'])
            if not hits: raise ValueError('no target alignment')
            stage='parse_m8CC'; rows=parse_m8(outputs['8CC'])
            stage='alignment_matching'; hit,row=select_primary(hits,rows,protein,dna)
            stage='traceback_validation'; vals,detail=verify(hit,row,protein,dna,r['region'])
            detail.update({'m10_alignment_count':len(hits),'m8CC_alignment_count':len(rows),'alignment_selection':'first reported m10 alignment; unique coordinate and traceback match'})
        r.update(vals); ev['tfasty']=detail
        ev['assessment_status']='tfasty_v3_full_reference' if detail['full_reference'] else 'tfasty_v3_partial_reference'
    except (KeyError,ValueError,TypeError,OSError,subprocess.SubprocessError) as exc:
        ev['assessment_status']='tfasty_v3_failed'
        ev['tfasty']={'error':str(exc),'failure_stage':stage,'error_type':type(exc).__name__}
        if 'outputs' in locals(): ev['tfasty']['raw_outputs']=outputs
        if getattr(exc,'stderr',None): ev['tfasty']['stderr']=str(exc.stderr)
    assess_boundaries(r, a.genomes, a.flank_window)
    vals=[r.get(k) for k in (*FIELDS,'size_mismatch')]
    r['known_disablements']=sum(v for v in vals if v is not None)
    r['total_disablements']=None if any(v is None for v in vals) else sum(vals)
    return r


def audit_row(r):
    e=r.get('mutation_evidence',{}); t=e.get('tfasty',{}); m=t.get('alignment_metadata',{})
    b=e.get('terminal_boundary',{}); bs=b.get('start',{}); bt=b.get('stop',{})
    def pct(k): return 100*float(m[k]) if k in m else None
    return dict(region_id=r['region']['region_id'],protein_id=r['best_protein']['protein_id'],
                tfasty_assessment_status=e.get('assessment_status'),tfasty_alignment_parsed='query_alignment' in t,
                tfasty_reference_coverage_pct=t.get('reference_coverage_pct'),tfasty_full_reference=t.get('full_reference'),
                tfasty_identity_pct=pct('sw_ident'),tfasty_similarity_pct=pct('sw_sim'),tfasty_bitscore=m.get('tfy_bits'),
                tfasty_evalue=m.get('tfy_expect'),tfasty_query_start=t.get('query_start'),tfasty_query_end=t.get('query_end'),
                tfasty_alignment_code=t.get('alignment_code'),m10_alignment_count=t.get('m10_alignment_count'),m8CC_alignment_count=t.get('m8CC_alignment_count'),failure_stage=t.get('failure_stage'),error=t.get('error'),
                boundary_status=b.get('status'),boundary_error=b.get('reason'),
                expected_start_status=bs.get('expected_site_status'),expected_stop_status=bt.get('expected_site_status'),
                expected_start_altered=bs.get('expected_site_altered'),expected_stop_altered=bt.get('expected_site_altered'),
                expected_start_codon=bs.get('expected_site',{}).get('codon'),expected_stop_codon=bt.get('expected_site',{}).get('codon'),
                expected_start_outside=bs.get('expected_site',{}).get('outside_candidate'),expected_stop_outside=bt.get('expected_site',{}).get('outside_candidate'),
                alternative_start_count=len(bs.get('alternatives',[])),alternative_stop_count=len(bt.get('alternatives',[])),
                **{k:r.get(k) for k in FIELDS})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path); p.add_argument('output',type=Path)
    p.add_argument('--genome', type=Path, required=True, help='original unmasked genomic FASTA')
    p.add_argument('--flank-window', type=int, default=30, help='alternative codon search radius in nt; multiple of 3')
    p.add_argument('--reference',type=Path,default=Path('uniprot_sprot.fasta'))
    p.add_argument('--tfasty',type=Path,default=Path('fasta-36.3.8i/bin/tfasty36'))
    p.add_argument('--limit',type=int); p.add_argument('--workers',type=int,default=1); p.add_argument('--timeout',type=int,default=120)
    a=p.parse_args()
    if a.input.resolve()==a.output.resolve() or a.workers<1 or a.timeout<1 or (a.limit is not None and a.limit<1): p.error('invalid arguments')
    if a.flank_window < 0 or a.flank_window % 3: p.error('--flank-window must be a nonnegative multiple of 3')
    if not a.genome.is_file(): p.error('genome FASTA unavailable')
    if not a.tfasty.is_file() or not os.access(a.tfasty,os.X_OK): p.error('TFASTY executable unavailable')
    records=[json.loads(l) for l in a.input.open() if l.strip()]; records=records[:a.limit]
    if not records:
        # An assembly with zero candidates is a valid workflow result.
        a.output.parent.mkdir(parents=True, exist_ok=True)
        paths = [a.output, a.output.with_suffix('.audit.tsv'),
                 a.output.with_suffix('.failed.jsonl')]
        if a.input.resolve() in [x.resolve() for x in paths]:
            p.error('output path conflicts with input')
        a.output.write_text('')
        paths[2].write_text('')
        dummy = {'region': {'region_id': ''}, 'best_protein': {'protein_id': ''}}
        with paths[1].open('w') as af:
            csv.DictWriter(af, fieldnames=list(audit_row(dummy)),
                           delimiter='\t').writeheader()
        print(json.dumps({'records': 0, 'parsed': 0}, indent=2))
        return
    if any(r.get(k) is not None for r in records for k in FIELDS): p.error('use original unassessed JSONL')
    a.genomes=load_refs(a.genome,{r['region']['chromosome'] for r in records})
    if not a.genomes: p.error('no candidate contig identifiers match genome FASTA')
    refs=load_refs(a.reference,{r['best_protein']['protein_id'] for r in records})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    audit=a.output.with_suffix('.audit.tsv'); failures=a.output.with_suffix('.failed.jsonl'); counts=Counter()
    paths=[a.output,audit,failures]; temps=[x.with_name(x.name+'.tmp') for x in paths]
    if a.input.resolve() in [x.resolve() for x in paths+temps]: p.error('output path conflicts with input')
    with ThreadPoolExecutor(max_workers=a.workers) as pool, temps[0].open('w') as out, temps[1].open('w') as af, temps[2].open('w') as ff:
        writer=csv.DictWriter(af,fieldnames=list(audit_row(records[0])),delimiter='\t'); writer.writeheader()
        for r in pool.map(lambda r:run_record(r,refs,a),records):
            counts['records']+=1; row=audit_row(r); writer.writerow(row); out.write(json.dumps(r)+'\n')
            counts['parsed' if row['tfasty_alignment_parsed'] else 'unassessed']+=1
            if not row['tfasty_alignment_parsed']: ff.write(json.dumps(r)+'\n'); counts['failure_'+str(row['failure_stage'])]+=1
            counts['boundary_'+str(row['boundary_status'])]+=1
            for kind in ('start','stop'):
                counts['expected_'+kind+'_assessed']+=row['expected_'+kind+'_altered'] is not None
                counts['expected_'+kind+'_altered']+=row['expected_'+kind+'_altered']==1
            for k in FIELDS:
                v=r.get(k); counts[k+'_assessed']+=v is not None; counts[k+'_positive']+=(v or 0)>0
    for tmp,path in zip(temps,paths): tmp.replace(path)
    print(json.dumps(counts,indent=2))

if __name__=='__main__': main()
