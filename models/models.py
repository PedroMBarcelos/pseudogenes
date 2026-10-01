from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass
class GenomicRegion:
    """Represents a merged genomic region."""

    region_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    parent_protein: str
    score_density: float
    num_hsps: int
    total_frameshifts: int
    internal_stops: int
    total_disablements: int

    @property
    def length(self) -> int:
        """Return region length."""
        return self.end - self.start + 1

    def to_dict(self) -> Dict[str, Any]:
        """Convert region to dictionary."""
        return {
            "region_id": self.region_id,
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "parent_protein": self.parent_protein,
            "score_density": self.score_density,
            "num_hsps": self.num_hsps,
            "total_frameshifts": self.total_frameshifts,
            "internal_stops": self.internal_stops,
            "total_disablements": self.total_disablements,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "GenomicRegion":
        """Create a GenomicRegion from a dictionary."""
        return cls(**d)


@dataclass
class BlastHit:
    """Represents an individual BLAST hit."""
    
    qseqid: str      # Query sequence ID
    qlen: int        # Query length
    sseqid: str      # Subject sequence ID
    slen: int        # Subject length
    qstart: int      # Query start position
    qend: int        # Query end position
    sstart: int      # Subject start position
    send: int        # Subject end position
    qseq: str        # Query aligned sequence
    sseq: str        # Subject aligned sequence
    evalue: float    # E-value
    bitscore: float  # Bit score
    score: int       # Raw score (BLOSUM62)
    length: int      # Alignment length
    pident: float    # Percent identity
    nident: int      # Number of identical matches
    mismatch: int    # Number of mismatches
    positive: int    # Number of positive matches
    gapopen: int     # Number of gap openings
    gaps: int        # Total gaps
    ppos: float      # Percent positive matches
    sframe: int      # Subject frame
    sstrand: str     # Subject strand
    qcovs: float     # Query coverage per subject
    qcovhsp: float   # Query coverage per HSP
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert hit to dictionary."""
        return {
            'qseqid': self.qseqid,
            'qlen': self.qlen,
            'sseqid': self.sseqid,
            'slen': self.slen,
            'qstart': self.qstart,
            'qend': self.qend,
            'sstart': self.sstart,
            'send': self.send,
            'qseq': self.qseq,
            'sseq': self.sseq,
            'evalue': self.evalue,
            'bitscore': self.bitscore,
            'score': self.score,
            'length': self.length,
            'pident': self.pident,
            'nident': self.nident,
            'mismatch': self.mismatch,
            'positive': self.positive,
            'gapopen': self.gapopen,
            'gaps': self.gaps,
            'ppos': self.ppos,
            'sframe': self.sframe,
            'sstrand': self.sstrand,
            'qcovs': self.qcovs,
            'qcovhsp': self.qcovhsp
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'BlastHit':
        return cls(**d)


@dataclass
class ProteinHitGroup:
    """Group of HSPs for a specific protein in a region."""
    
    protein_id: str           # ID da proteína (qseqid)
    hsps: List[BlastHit]      # Lista de HSPs
    total_score: int          # Soma dos scores
    total_length: int         # Soma dos comprimentos das sequências
    score_density: float      # Densidade de score
    
    @classmethod
    def from_hsps(cls, protein_id: str, hsps: List[BlastHit]) -> 'ProteinHitGroup':
        """
        Create a ProteinHitGroup from a list of HSPs.
        
        Follows expert specification:
        Score Density = Σ(Scores of HSPs) / Σ(Alignment Lengths of HSPs)
        
        Where:
        - Numerator: Sum of raw scores from each HSP
        - Denominator: Sum of actual alignment lengths (hsp.length)
          NOT query coverage (qend-qstart) which would include gaps between HSPs
        
        Args:
            protein_id: Protein ID
            hsps: List of HSPs for this protein
        
        Returns:
            ProteinHitGroup object with calculated density
        """
        total_score = sum(hsp.score for hsp in hsps)
        # Use actual alignment length (number of aligned residues)
        # NOT query coverage which includes gaps between HSPs
        total_length = sum(hsp.length for hsp in hsps)
        
        # Calculate density: score / alignment_length
        score_density = total_score / total_length if total_length > 0 else 0.0
        
        return cls(
            protein_id=protein_id,
            hsps=hsps,
            total_score=total_score,
            total_length=total_length,
            score_density=score_density
        )
    
    def to_dict(self) -> Dict:
        return {
            'protein_id': self.protein_id,
            'hsps': [hsp.to_dict() for hsp in self.hsps],
            'total_score': self.total_score,
            'total_length': self.total_length,
            'score_density': self.score_density,
        }

    @classmethod
    def from_dict(cls, d: Dict) -> 'ProteinHitGroup':
        return cls(
            protein_id=d['protein_id'],
            hsps=[BlastHit.from_dict(h) for h in d['hsps']],
            total_score=d['total_score'],
            total_length=d['total_length'],
            score_density=d['score_density'],
        )


@dataclass
class RegionAnnotation:
    """Complete annotation of a genomic region."""
    
    region: GenomicRegion           # Região genômica fundida
    best_protein: ProteinHitGroup   # Proteína com maior densidade de score
    all_proteins: List[ProteinHitGroup]  # Todas as proteínas que mapearam

    def to_json_dict(self) -> Dict:
        """Serialização completa para JSON."""
        return {
            'region': self.region.to_dict(),
            'best_protein': self.best_protein.to_dict(),
            'all_proteins': [p.to_dict() for p in self.all_proteins],
        }
    
    def to_dict(self) -> Dict:
        """Convert annotation to dictionary for TSV files."""
        return {
            'chromosome': self.region.chromosome,
            'start': self.region.start,
            'end': self.region.end,
            'length': self.region.length,
            'strand': self.region.strand,
            'best_protein_id': self.best_protein.protein_id,
            'best_protein_score_density': self.best_protein.score_density,
            'best_protein_total_score': self.best_protein.total_score,
            'best_protein_total_length': self.best_protein.total_length,
            'best_protein_num_hsps': len(self.best_protein.hsps),
            'num_competing_proteins': len(self.all_proteins)
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'RegionAnnotation':
        """Deserializa um RegionAnnotation a partir de um dicionário."""
        return cls(
            region=GenomicRegion.from_dict(d['region']),
            best_protein=ProteinHitGroup.from_dict(d['best_protein']),
            all_proteins=[ProteinHitGroup.from_dict(p) for p in d['all_proteins']],
        )


@dataclass
class DisablementCounts:
    """Contadores de diferentes tipos de mutações que podem inativar um gene."""
    
    non_synonymous_substitutions: int = 0  # Substituições não-sinônimas
    in_frame_indels: int = 0               # Indels in-frame
    frameshifts: int = 0                   # Mudanças de frame
    missing_start_codon: int = 0           # Perda do start codon
    missing_stop_codon: int = 0            # Perda do stop codon
    premature_stop_codons: int = 0         # Stop codons prematuros
    size_mismatch: int = 0                 # Tamanho incompatível (Regra dos 20%)
    
    @property
    def total_disablements(self) -> int:
        """Retorna o total de mutações inativadoras."""
        return (
            self.non_synonymous_substitutions +
            self.in_frame_indels +
            self.frameshifts +
            self.missing_start_codon +
            self.missing_stop_codon +
            self.premature_stop_codons +
            self.size_mismatch
        )
    
    def to_dict(self) -> Dict[str, int]:
        """Converte os contadores em dicionário."""
        return {
            "non_synonymous_substitutions": self.non_synonymous_substitutions,
            "in_frame_indels": self.in_frame_indels,
            "frameshifts": self.frameshifts,
            "missing_start_codon": self.missing_start_codon,
            "missing_stop_codon": self.missing_stop_codon,
            "premature_stop_codons": self.premature_stop_codons,
            "size_mismatch": self.size_mismatch,
            "total_disablements": self.total_disablements
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'DisablementCounts':
        return cls(
            non_synonymous_substitutions=d['non_synonymous_substitutions'],
            in_frame_indels=d['in_frame_indels'],
            frameshifts=d['frameshifts'],
            missing_start_codon=d['missing_start_codon'],
            missing_stop_codon=d['missing_stop_codon'],
            premature_stop_codons=d['premature_stop_codons'],
            size_mismatch=d['size_mismatch'],
            # total_disablements ignorado — é @property calculada automaticamente
        )


@dataclass
class PseudogeneAnnotation:
    """
    Anotação completa de um pseudogene potencial.
    
    Estrutura de acesso aos atributos:
    - annotation.is_pseudogene
    - annotation.region_annotation.region.chromosome
    - annotation.region_annotation.region.start
    - annotation.region_annotation.region.end
    - annotation.region_annotation.region.strand
    - annotation.region_annotation.best_protein.protein_id
    - annotation.region_annotation.best_protein.score_density
    - annotation.disablements.non_synonymous_substitutions
    - annotation.disablements.in_frame_indels
    - annotation.disablements.frameshifts
    - annotation.disablements.premature_stop_codons
    - annotation.disablements.missing_start_codon
    - annotation.disablements.missing_stop_codon
    - annotation.disablements.total_disablements (propriedade)
    """
    
    region_annotation: RegionAnnotation  # Anotação da região
    disablements: DisablementCounts      # Contadores de mutações
    is_pseudogene: bool                  # Se é classificado como pseudogene
    is_small_orf: bool = False           # If region is small (<300 nt / 100 aa)
    
    def to_json_dict(self):
        """Converte a anotação em dicionário."""
        result = self.region_annotation.to_json_dict()
        result.update(self.disablements.to_dict())
        result['is_pseudogene'] = self.is_pseudogene
        result['is_small_orf'] = self.is_small_orf
        
        # Add concatenated sequences for visualization
        hsps = self.region_annotation.best_protein.hsps
        if hsps:
            # Sort HSPs by query start to ensure correct order
            sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
            result['qseq'] = "".join(h.qseq for h in sorted_hsps)
            result['sseq'] = "".join(h.sseq for h in sorted_hsps)
        else:
            result['qseq'] = ""
            result['sseq'] = ""
            
        return result
    
    def to_tsv_dict(self):
        result = self.region_annotation.to_dict()
        result.update(self.disablements.to_dict())
        result['is_pseudogene'] = self.is_pseudogene
        result['is_small_orf'] = self.is_small_orf
        
        # Add concatenated sequences for visualization
        hsps = self.region_annotation.best_protein.hsps
        if hsps:
            # Sort HSPs by query start to ensure correct order
            sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
            result['qseq'] = "".join(h.qseq for h in sorted_hsps)
            result['sseq'] = "".join(h.sseq for h in sorted_hsps)
        else:
            result['qseq'] = ""
            result['sseq'] = ""
            
        return result

    @classmethod
    def from_dict(cls, d: Dict) -> 'PseudogeneAnnotation':
        # Reconstrói RegionAnnotation a partir da estrutura plana
        region = GenomicRegion(
            region_id=d["region"]["region_id"],
            chromosome=d["region"]["chromosome"],
            start=d["region"]["start"],
            end=d["region"]["end"],
            strand=d["region"]["strand"],
            parent_protein=d["region"]["parent_protein"],
            score_density=d["region"]["score_density"],
            num_hsps=d["region"]["num_hsps"],
            total_frameshifts=d["region"]["total_frameshifts"],
            internal_stops=d["region"]["internal_stops"],
            total_disablements=d["region"]["total_disablements"],
        )
        region_annotation = RegionAnnotation(
            region=region,
            best_protein=ProteinHitGroup.from_dict(d['best_protein']),
            all_proteins=[ProteinHitGroup.from_dict(p) for p in d['all_proteins']],
        )
        disablements = DisablementCounts(
            non_synonymous_substitutions=d['non_synonymous_substitutions'],
            in_frame_indels=d['in_frame_indels'],
            frameshifts=d['frameshifts'],
            missing_start_codon=d['missing_start_codon'],
            missing_stop_codon=d['missing_stop_codon'],
            premature_stop_codons=d['premature_stop_codons'],
            size_mismatch=d['size_mismatch'],
        )
        return cls(
            region_annotation=region_annotation,
            disablements=disablements,
            is_pseudogene=d['is_pseudogene'],
            is_small_orf=d.get('is_small_orf', False),
        )