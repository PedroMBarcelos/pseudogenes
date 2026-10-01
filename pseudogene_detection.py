import json
import logging
from pathlib import Path
from typing import Dict, List, Set
from Bio import SeqIO


from models.models import BlastHit, ProteinHitGroup, PseudogeneAnnotation, RegionAnnotation, DisablementCounts


logger = logging.getLogger(__name__)

# Caractere que representa gap em alinhamentos
GAP_CHAR: str = "-"

# Caractere que representa stop codon em sequências traduzidas
STOP_CODON_CHAR: str = "*"

# Aminoácido inicial esperado (Metionina)
START_CODON_AA: str = "M"

# Limites de tamanho para a Regra dos 20%
MIN_LENGTH_RATIO: float = 0.8  # 80% do tamanho da referência
MAX_LENGTH_RATIO: float = 1.2  # 120% do tamanho da referência

# Minimum region size for Small ORF classification (nucleotides)
MIN_REGION_SIZE_NT: int = 300  # Regions below this are classified as 'small ORF' (100 aa)
MIN_ALIGNMENT_SIZE_AA: int = 100  # Alignment length threshold for Small ORF check


def calculate_subject_coverage_nt(hsps: List[BlastHit]) -> int:
    """
    Calcula a cobertura no SUBJECT (Genoma) em nucleotídeos, excluindo gaps.
    
    Args:
        hsps: Lista de HSPs
        
    Returns:
        Total de nucleotídeos cobertos pelo alinhamento
    """
    if not hsps:
        return 0
    
    # Coletar intervalos no genoma (Subject)
    # Nota: tblastn pode ter sstart > send se for na fita reversa
    intervals = []
    for hsp in hsps:
        start = min(hsp.sstart, hsp.send)
        end = max(hsp.sstart, hsp.send)
        intervals.append((start, end))
    
    # Ordenar por posição inicial
    intervals.sort()
    
    # Fundir intervalos sobrepostos
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    
    # Somar comprimentos
    total = sum(end - start + 1 for start, end in merged)
    
    return total


def count_non_synonymous_substitutions(query_seq: str, subject_seq: str) -> int:
    """
    Conta substituições não-sinônimas entre sequências alinhadas.
    
    Conta aminoácidos não idênticos na mesma posição do alinhamento,
    excluindo gaps. IMPORTANTE: Nem toda substituição indica pseudogene,
    apenas contamos para análise.
    
    Args:
        query_seq: Sequência query alinhada (proteína funcional)
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de substituições não-sinônimas
    """
    if len(query_seq) != len(subject_seq):
        logger.warning(
            f"Sequências com comprimentos diferentes: "
            f"query={len(query_seq)}, subject={len(subject_seq)}"
        )
    
    count = 0
    min_length = min(len(query_seq), len(subject_seq))
    
    for i in range(min_length):
        q_aa = query_seq[i]
        s_aa = subject_seq[i]
        
        # Ignora posições com gaps
        if q_aa == GAP_CHAR or s_aa == GAP_CHAR:
            continue
        
        # Ignora stop codons na contagem de substituições
        if q_aa == STOP_CODON_CHAR or s_aa == STOP_CODON_CHAR:
            continue
        
        # Conta se forem diferentes
        if q_aa != s_aa:
            count += 1
    
    # Calcula percentual de diferença
    total_positions = min_length - query_seq.count(GAP_CHAR) - subject_seq.count(GAP_CHAR)
    if total_positions > 0:
        percent_diff = (count / total_positions) * 100
        logger.debug(f"Substituições: {count}/{total_positions} ({percent_diff:.1f}%)")
    
    return count


def count_in_frame_indels(query_seq: str, subject_seq: str) -> int:
    """
    Conta indels in-frame (gaps individuais) presentes em ambas as sequências.
    
    Args:
        query_seq: Sequência query alinhada
        subject_seq: Sequência subject alinhada
    
    Returns:
        Número total de gaps individuais
    """
    query_gaps = query_seq.count(GAP_CHAR)
    subject_gaps = subject_seq.count(GAP_CHAR)
    
    total_gaps = query_gaps + subject_gaps
    
    return total_gaps


def count_frameshifts(hsps: List[BlastHit]) -> int:
    """
    Calcula o número de frameshifts baseado em frames distintos dos HSPs.
    
    Frameshifts = count(distinct_frames) - 1
    
    Args:
        hsps: Lista de HSPs da mesma sequência de referência
    
    Returns:
        Número de frameshifts
    """
    if not hsps:
        return 0
    
    # Coleta todos os frames únicos
    frames: Set[int] = set()
    for hsp in hsps:
        frames.add(hsp.sframe)
    
    # Número de frameshifts é o número de frames distintos menos 1
    num_frameshifts = len(frames) - 1
    
    logger.debug(
        f"Frames detectados: {sorted(frames)}, "
        f"Frameshifts: {num_frameshifts}"
    )
    
    return max(0, num_frameshifts)


def check_missing_start_codon(subject_seq: str) -> bool:
    """
    Verifica se o start codon (Metionina) está ausente.
    
    IMPORTANTE: tblastn traduz o DNA mas NÃO inclui o start codon (ATG->M) 
    automaticamente na sequência traduzida. Por isso, esta verificação é
    relevante apenas se o alinhamento começa no início real do gene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o start codon está ausente, False caso contrário
    """
    # Remove gaps do início
    trimmed_seq = subject_seq.lstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, a ausência de M no início pode ser normal se o alinhamento
    # não começa exatamente no start codon. Portanto, NÃO consideramos isso
    # como evidência de pseudogene por padrão.
    # Apenas reportamos se estiver completamente ausente E tiver outros indicadores
    first_aa = trimmed_seq[0]
    
    # Relaxa: só marca como problemático se não for M E a sequência for longa o suficiente
    # para esperar que comece no início real
    if len(trimmed_seq) < 50:  # Alinhamentos curtos podem não incluir o start
        return False
    
    is_missing = first_aa != START_CODON_AA
    
    if is_missing:
        logger.debug(f"Possível start codon ausente. Primeiro AA: {first_aa}")
    
    # NÃO marca como problema por padrão - retorna False
    return False


def check_missing_stop_codon(subject_seq: str) -> bool:
    """
    Verifica se o stop codon está ausente no final da sequência.
    
    IMPORTANTE: tblastn traduz até encontrar um stop codon, mas NÃO inclui
    o stop codon (*) na sequência traduzida resultante. Portanto, a ausência
    de * no final é NORMAL e NÃO indica pseudogene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o stop codon está ausente, False caso contrário
    """
    # Remove gaps do final
    trimmed_seq = subject_seq.rstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, o stop codon NÃO aparece na sequência traduzida normalmente.
    # A presença de * no final seria anormal. A ausência é esperada.
    # Portanto, SEMPRE retornamos False (não é problema)
    last_aa = trimmed_seq[-1]
    
    # Se houver * no final, isso seria estranho, mas não necessariamente problema
    if last_aa == STOP_CODON_CHAR:
        logger.debug(f"Stop codon presente no final (incomum mas não problemático): {last_aa}")
    
    # NÃO marca como problema - retorna sempre False
    return False


def count_premature_stop_codons(subject_seq: str) -> int:
    """
    Conta stop codons prematuros (internos) na sequência.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de stop codons internos
    """
    # Remove gaps
    clean_seq = subject_seq.replace(GAP_CHAR, '')
    
    if not clean_seq:
        return 0
    
    # Remove o último caractere (que pode ser stop codon legítimo)
    internal_seq = clean_seq[:-1] if len(clean_seq) > 1 else ""
    
    # Conta stop codons internos
    count = internal_seq.count(STOP_CODON_CHAR)
    
    if count > 0:
        logger.debug(f"Stop codons prematuros detectados: {count}")
    
    return count


def analyze_protein_for_disablements(
    protein_hit_group: ProteinHitGroup
) -> DisablementCounts:
    """
    Analyze HSP group for inactivating mutations.
    
    Args:
        protein_hit_group: Protein HSP group
    
    Returns:
        Counters for different mutation types
    """
    logger.debug(f"Analyzing protein {protein_hit_group.protein_id}")
    
    counts = DisablementCounts()
    
    hsps = protein_hit_group.hsps
    
    if not hsps:
        return counts
    
    # Concatenate all HSP sequences for global analysis
    all_query_seq = ""
    all_subject_seq = ""
    
    for hsp in hsps:
        all_query_seq += hsp.qseq
        all_subject_seq += hsp.sseq
    
    # 1. Non-synonymous substitutions
    counts.non_synonymous_substitutions = count_non_synonymous_substitutions(
        all_query_seq, all_subject_seq
    )
    
    # 2. In-frame indels
    counts.in_frame_indels = count_in_frame_indels(
        all_query_seq, all_subject_seq
    )
    
    # 3. Frameshifts
    counts.frameshifts = count_frameshifts(hsps)
    
    # 4. Missing start codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_start_codon = 1 if check_missing_start_codon(all_subject_seq) else 0
    
    # 5. Missing stop codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_stop_codon = 1 if check_missing_stop_codon(all_subject_seq) else 0
    
    # 6. Premature stop codons
    counts.premature_stop_codons = count_premature_stop_codons(all_subject_seq)
    
    
    return counts


def check_start_stop_codons_genomic(
    genome_seqs: Dict,  # Accepts the preloaded dictionary directly
    chromosome: str,
    start: int,
    end: int,
    strand: str,
    reference_length_aa: int,
    check_start: bool = True,
    check_stop: bool = True
) -> Dict[str, bool]:
    """
    Check for start/stop codons by expanding aligned region in preloaded genome dictionary.
    """
    if chromosome not in genome_seqs:
        logger.warning(f"Chromosome {chromosome} not found in genome")
        return {"has_start": False, "has_stop": False, "start_codon": None, "stop_codon": None}
    
    genome_seq = genome_seqs[chromosome]
    
    # Calculate dynamic padding (max 20% of reference length)
    max_padding_nt = int((reference_length_aa * 3) * 0.20)
    padding = max_padding_nt
    
    # Convert to 0-based
    start_0 = start - 1
    end_0 = end
    
    # Adjust coordinates with padding
    search_start = max(0, start_0 - padding)
    search_end = min(len(genome_seq), end_0 + padding)
    
    # Extract extended region
    raw_seq = genome_seq[search_start:search_end]
    
    # Orient sequence (reverse complement if negative strand)
    if strand == '-':
        raw_seq = raw_seq.reverse_complement()
    
    dna_str = str(raw_seq).upper()
    
    # Calculate relative alignment position within extracted string
    if strand == '+':
        rel_align_start = start_0 - search_start
        rel_align_end = end_0 - search_start
    else:
        rel_align_start = search_end - end_0
        rel_align_end = search_end - start_0
    
    valid_starts = ["ATG", "GTG", "TTG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # --- START CODON SEARCH (UPSTREAM) ---
    has_start = False
    start_codon = None
    
    if check_start:
        upstream_region = dna_str[0:rel_align_start]
        for i in range(len(upstream_region) - 3, -1, -3):
            codon = upstream_region[i:i+3]
            if len(codon) < 3:
                continue
            if codon in valid_starts:
                has_start = True
                start_codon = codon
                break
            if codon in stop_codons:
                break
    
    # --- STOP CODON SEARCH (DOWNSTREAM) ---
    has_stop = False
    stop_codon = None
    
    if check_stop:
        downstream_region = dna_str[rel_align_end:]
        for i in range(0, len(downstream_region) - 2, 3):
            codon = downstream_region[i:i+3]
            if len(codon) < 3:
                continue
            if codon in stop_codons:
                has_stop = True
                stop_codon = codon
                break
    
    return {
        "has_start": has_start,
        "has_stop": has_stop,
        "start_codon": start_codon,
        "stop_codon": stop_codon
    }

def classify_as_pseudogene(
    disablements: DisablementCounts,
    min_disablements: int = 1
) -> bool:
    """
    Classifica se uma região é um pseudogene baseado no número de mutações.
    
    CRITÉRIOS ATUALIZADOS:
    - Stop codons prematuros são forte evidência de pseudogene
    - Frameshifts são forte evidência de pseudogene
    - Ausência de start codon (verificada genomicamente) é forte evidência
    - Ausência de stop codon (verificada genomicamente) é forte evidência
    - Tamanho <80% ou >120% da proteína de referência é forte evidência (Regra dos 20%)
    - Substituições e indels sozinhos NÃO indicam pseudogene (podem ser variação normal)
    
    Args:
        disablements: Contadores de mutações
        min_disablements: Número mínimo de mutações para classificar como pseudogene
    
    Returns:
        True se for classificado como pseudogene, False caso contrário
    """
    # Evidências fortes de pseudogene
    strong_evidence = (
        disablements.premature_stop_codons +
        disablements.frameshifts +
        disablements.missing_start_codon +
        disablements.missing_stop_codon +
        disablements.size_mismatch
    )
    
    # Classifica como pseudogene apenas se houver evidência forte
    is_pseudogene = strong_evidence >= min_disablements
    
    logger.debug(
        f"Classificação: strong_evidence={strong_evidence}, "
        f"is_pseudogene={is_pseudogene}"
    )
    
    return is_pseudogene


def annotate_pseudogenes(
    region_annotations: List[RegionAnnotation],
    genome_fasta_path: Path,
    output_path: Path,
    min_disablements: int = 1,
) -> List[PseudogeneAnnotation]:
    
    # FIX 1: Load genome ONCE outside the loop to prevent massive I/O bottlenecks
    logger.info("Loading genome FASTA into memory...")
    genome_seqs = {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta_path, "fasta")}
    
    pseudogene_annotations = []
    
    for region_ann in region_annotations:
        disablements = analyze_protein_for_disablements(region_ann.best_protein)

        hsps = region_ann.best_protein.hsps
        sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
        first_hsp = sorted_hsps[0]
        last_hsp = sorted_hsps[-1]
        
        first_qseq_clean = first_hsp.qseq.lstrip('-')
        has_start_in_alignment = (
            len(first_qseq_clean) > 0 and
            first_qseq_clean[0] == 'M'
        )
        has_stop_in_alignment = (last_hsp.qend == last_hsp.qlen)
        
        need_genomic_check = not has_start_in_alignment or not has_stop_in_alignment
        
        if need_genomic_check:
            genomic_check = check_start_stop_codons_genomic(
                genome_seqs,  
                region_ann.region.chromosome,
                region_ann.region.start,
                region_ann.region.end,
                region_ann.region.strand,
                reference_length_aa=region_ann.best_protein.hsps[0].qlen,
                check_start=not has_start_in_alignment,
                check_stop=not has_stop_in_alignment
            )
            
            disablements.missing_start_codon = 0 if genomic_check["has_start"] else 1
            disablements.missing_stop_codon = 0 if genomic_check["has_stop"] else 1
        else:
            disablements.missing_start_codon = 0
            disablements.missing_stop_codon = 0

        # FIX 2: Implement the 20% size rule instead of hardcoding 0
        ref_len_aa = first_hsp.qlen
        alignment_len_aa = region_ann.best_protein.total_length # Or appropriate aa length metric
        if ref_len_aa > 0:
            ratio = alignment_len_aa / ref_len_aa
            if ratio < MIN_LENGTH_RATIO or ratio > MAX_LENGTH_RATIO:
                disablements.size_mismatch = 1
            else:
                disablements.size_mismatch = 0
        else:
            disablements.size_mismatch = 0
        
        is_normal_size = (region_ann.region.length >= MIN_REGION_SIZE_NT) or (alignment_len_aa >= MIN_ALIGNMENT_SIZE_AA)
        is_small_orf = not is_normal_size
        
        is_pseudogene = classify_as_pseudogene(disablements, min_disablements)
        
        pseudogene_annotations.append(
            PseudogeneAnnotation(
                region_annotation=region_ann,
                disablements=disablements,
                is_pseudogene=is_pseudogene,
                is_small_orf=is_small_orf
            )
        )
        

    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.debug(f"Pseudogenes: {num_pseudogenes}/{len(pseudogene_annotations)}")

    logger.debug(f"Saving {len(pseudogene_annotations)} pseudogene annotations to: {output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in pseudogene_annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha

    return pseudogene_annotations


def load_region_annotations_from_json(annotations_path: Path) -> List[RegionAnnotation]:
    """
    Carrega RegionAnnotations de um arquivo JSON Lines.

    Args:
        annotations_path: Caminho para o arquivo .jsonl gerado por annotate_regions_with_best_proteins

    Returns:
        Lista de RegionAnnotation completamente reconstruída
    """
    annotations = []
    with annotations_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line:
                d = json.loads(line)
                annotations.append(RegionAnnotation.from_dict(d))

    logger.debug(f"Loaded {len(annotations)} region annotations from {annotations_path}")
    return annotations

import sys
if __name__ == "__main__":
    logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
        )
    region_annotations = load_region_annotations_from_json(Path(sys.argv[1]))
    annotate_pseudogenes(region_annotations, Path(sys.argv[2]), Path(sys.argv[3]))