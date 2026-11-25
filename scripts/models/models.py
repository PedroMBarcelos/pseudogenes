
from dataclasses import dataclass
from typing import Any, Dict, List, Optional


# ============================================================
# GENOMIC REGION
# ============================================================

@dataclass
class GenomicRegion:
    """
    Represents a merged genomic region.

    Coordinates are assumed to follow the BED convention:
        - 0-based start
        - end-exclusive
    """

    region_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    parent_protein: str
    score_density: float
    num_hsps: int
    total_frameshifts: Optional[int]
    internal_stops: Optional[int]
    total_disablements: Optional[int]

    @property
    def length(self) -> int:
        """Return genomic region length in nucleotides."""
        return self.end - self.start

    def to_dict(self) -> Dict[str, Any]:
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
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "GenomicRegion":
        return cls(**d)


# ============================================================
# TBLASTN HIT
# ============================================================

@dataclass
class BlastHit:
    """
    Represents an original TBLASTN HSP.

    IMPORTANT:
    These values come from TBLASTN, not SSEARCH.

    They are retained for genomic mapping and descriptive
    analysis of the original translated alignments.
    """

    qseqid: str
    qlen: int
    sseqid: str
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    qseq: str
    sseq: str
    evalue: float
    bitscore: float
    score: int
    length: int
    pident: float
    nident: int
    mismatch: int
    positive: int
    gapopen: int
    gaps: int
    ppos: float
    sframe: int
    sstrand: str
    qcovs: float
    qcovhsp: float

    def to_dict(self) -> Dict[str, Any]:
        return vars(self).copy()

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "BlastHit":
        return cls(**d)


# ============================================================
# PROTEIN HIT GROUP
# ============================================================

@dataclass
class ProteinHitGroup:
    """
    Group of TBLASTN HSPs associated with one protein.

    The score values in this class are TBLASTN-derived.

    SSEARCH statistics are preserved separately inside
    RegionAnnotation.ssearch_selection.

    Score density:
        sum(raw TBLASTN scores) / sum(HSP alignment lengths)

    Note:
        Overlapping HSPs may contribute repeated evidence.
        This metric is descriptive and should not replace
        a properly calculated alignment coverage metric.
    """

    protein_id: str
    hsps: List[BlastHit]
    total_score: int
    total_length: int
    score_density: float

    @classmethod
    def from_hsps(
        cls,
        protein_id: str,
        hsps: List[BlastHit]
    ) -> "ProteinHitGroup":

        total_score = sum(
            hit.score for hit in hsps
        )

        total_length = sum(
            hit.length for hit in hsps
        )

        score_density = (
            total_score / total_length
            if total_length > 0
            else 0.0
        )

        return cls(
            protein_id=protein_id,
            hsps=hsps,
            total_score=total_score,
            total_length=total_length,
            score_density=score_density,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "protein_id": self.protein_id,
            "hsps": [
                hit.to_dict()
                for hit in self.hsps
            ],
            "total_score": self.total_score,
            "total_length": self.total_length,
            "score_density": self.score_density,
        }

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "ProteinHitGroup":

        return cls(
            protein_id=d["protein_id"],
            hsps=[
                BlastHit.from_dict(hit)
                for hit in d["hsps"]
            ],
            total_score=d["total_score"],
            total_length=d["total_length"],
            score_density=d["score_density"],
        )


# ============================================================
# REGION ANNOTATION
# ============================================================

@dataclass
class RegionAnnotation:
    """
    Complete annotation of a candidate genomic region.

    best_protein:
        Parent protein selected using SSEARCH ranking.

    best_protein.hsps:
        Original TBLASTN alignments belonging to that protein.

    ssearch_selection:
        Actual SSEARCH statistics used for parent selection.

    These data sources must remain distinguishable.
    """

    region: GenomicRegion
    best_protein: ProteinHitGroup
    all_proteins: List[ProteinHitGroup]

    ssearch_selection: Optional[Dict[str, Any]] = None

    def to_json_dict(self) -> Dict[str, Any]:

        result = {
            "region": self.region.to_dict(),
            "best_protein": self.best_protein.to_dict(),
            "all_proteins": [
                protein.to_dict()
                for protein in self.all_proteins
            ],
        }

        if self.ssearch_selection is not None:
            result["ssearch_selection"] = (
                self.ssearch_selection
            )

        return result

    def to_dict(self) -> Dict[str, Any]:
        """
        Flat representation for TSV output.

        Legacy best_protein score fields remain
        TBLASTN-derived for compatibility.
        """

        result = {
            "region_id": self.region.region_id,
            "chromosome": self.region.chromosome,
            "start": self.region.start,
            "end": self.region.end,
            "length": self.region.length,
            "strand": self.region.strand,

            "best_protein_id":
                self.best_protein.protein_id,

            "best_protein_score_density":
                self.best_protein.score_density,

            "best_protein_total_score":
                self.best_protein.total_score,

            "best_protein_total_length":
                self.best_protein.total_length,

            "best_protein_num_hsps":
                len(self.best_protein.hsps),

            "num_competing_proteins":
                len(self.all_proteins),
        }

        if self.ssearch_selection is not None:

            selection = self.ssearch_selection

            result.update({
                "ssearch_score_density":
                    selection.get("score_density"),

                "ssearch_total_bitscore":
                    selection.get("total_bitscore"),

                "ssearch_aligned_length":
                    selection.get("aligned_length"),

                "ssearch_hsp_count":
                    selection.get("hsp_count"),
            })

        return result

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "RegionAnnotation":

        annotation = cls(
            region=GenomicRegion.from_dict(
                d["region"]
            ),

            best_protein=ProteinHitGroup.from_dict(
                d["best_protein"]
            ),

            all_proteins=[
                ProteinHitGroup.from_dict(p)
                for p in d["all_proteins"]
            ],

            ssearch_selection=d.get(
                "ssearch_selection"
            ),
        )

        # Verify consistency of SSEARCH selection.
        if annotation.ssearch_selection is not None:

            selected_id = (
                annotation.ssearch_selection.get(
                    "best_protein"
                )
            )

            if (
                selected_id is not None
                and selected_id
                != annotation.best_protein.protein_id
            ):
                raise ValueError(
                    f"SSEARCH selected {selected_id}, "
                    f"but best_protein is "
                    f"{annotation.best_protein.protein_id} "
                    f"for region "
                    f"{annotation.region.region_id}"
                )

        return annotation


# ============================================================
# DISABLEMENT COUNTS
# ============================================================

@dataclass
class DisablementCounts:
    """
    Tracks descriptive alignment differences and
    potentially gene-inactivating mutations.

    IMPORTANT:
        None means the measurement was not assessed.

        Zero means the measurement was performed and
        no event was detected.

    Ordinary amino-acid substitutions and protein
    alignment gaps are NOT considered disablements
    by themselves.
    """

    # Descriptive protein-alignment differences.
    # Legacy field names retained for compatibility.
    non_synonymous_substitutions: int = 0
    in_frame_indels: int = 0

    # Potential gene-disrupting events.
    # Only populate with independently validated evidence.
    frameshifts: Optional[int] = None
    missing_start_codon: Optional[int] = None
    missing_stop_codon: Optional[int] = None
    premature_stop_codons: Optional[int] = None
    size_mismatch: Optional[int] = None

    def _disruption_values(self):
        return (
            self.frameshifts,
            self.missing_start_codon,
            self.missing_stop_codon,
            self.premature_stop_codons,
            self.size_mismatch,
        )

    @property
    def known_disablements(self) -> int:
        """
        Sum only assessed disruption counts.

        This is NOT necessarily the complete number
        of disablements.
        """

        return sum(
            value
            for value in self._disruption_values()
            if value is not None
        )

    @property
    def total_disablements(self) -> Optional[int]:
        """
        Total number of validated disablements.

        Returns None if any required disruption
        assessment is still missing.
        """

        values = self._disruption_values()

        if any(value is None for value in values):
            return None

        return sum(values)

    def to_dict(self) -> Dict[str, Any]:

        return {
            "non_synonymous_substitutions":
                self.non_synonymous_substitutions,

            "in_frame_indels":
                self.in_frame_indels,

            "frameshifts":
                self.frameshifts,

            "missing_start_codon":
                self.missing_start_codon,

            "missing_stop_codon":
                self.missing_stop_codon,

            "premature_stop_codons":
                self.premature_stop_codons,

            "size_mismatch":
                self.size_mismatch,

            "known_disablements":
                self.known_disablements,

            "total_disablements":
                self.total_disablements,
        }

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "DisablementCounts":

        return cls(
            non_synonymous_substitutions=d.get(
                "non_synonymous_substitutions", 0
            ),

            in_frame_indels=d.get(
                "in_frame_indels", 0
            ),

            frameshifts=d.get(
                "frameshifts"
            ),

            missing_start_codon=d.get(
                "missing_start_codon"
            ),

            missing_stop_codon=d.get(
                "missing_stop_codon"
            ),

            premature_stop_codons=d.get(
                "premature_stop_codons"
            ),

            size_mismatch=d.get(
                "size_mismatch"
            ),
        )


# ============================================================
# CANDIDATE REGION ANNOTATION
# ============================================================

@dataclass
class PseudogeneAnnotation:
    """
    Protein-guided genomic candidate characterization.

    Classification has only two possible values:

        small_orf:
            A short complete ORF supported by the
            selected parent-protein alignment.

        candidate:
            A protein-associated genomic region that
            does not meet the small-ORF criteria.

    Classification is independent of:
        - Reference protein coverage
        - Number of substitutions
        - Number of possible frameshifts
        - Number of possible internal stops

    None means a measurement was not assessed.

    The historical name PseudogeneAnnotation is
    retained for compatibility with existing scripts.
    """

    # Existing fields
    region_annotation: RegionAnnotation
    disablements: DisablementCounts

    # Classification
    classification_status: str = "candidate"
    is_small_orf: bool = False

    # Legacy field: no confirmed pseudogene classification
    is_pseudogene: Optional[bool] = None

    # Genomic sequence characterization
    # Sequence is oriented according to region.strand.
    sequence_nt: Optional[str] = None
    region_length_nt: Optional[int] = None

    # Coverage is descriptive, never a filter.
    reference_coverage: Optional[float] = None

    # A short region is not automatically a small ORF.
    is_small_region: Optional[bool] = None

    # Complete ORFs identified in this genomic sequence.
    # Each entry may contain:
    # start, end, frame, length_nt, sequence_nt,
    # and protein_support.
    orfs: Optional[List[Dict[str, Any]]] = None

    best_orf_length_nt: Optional[int] = None

    # Possible mutation observations.
    # These do not determine the classification.
    mutation_evidence: Optional[Dict[str, Any]] = None

    # --------------------------------------------------------
    # DERIVED PROPERTIES
    # --------------------------------------------------------

    @property
    def best_protein_id(self) -> str:
        return self.region_annotation.best_protein.protein_id

    @property
    def num_orfs(self) -> int:
        return len(self.orfs or [])

    # --------------------------------------------------------
    # JSON SERIALIZATION
    # --------------------------------------------------------

    def to_json_dict(self) -> Dict[str, Any]:
        """
        Serialize the candidate into the existing JSONL
        structure, including the new characterization fields.
        """

        result = self.region_annotation.to_json_dict()

        # Retain the descriptive and disruption measurements.
        result.update(self.disablements.to_dict())

        result.update({
            "classification_status":
                self.classification_status,

            "is_small_orf":
                self.is_small_orf,

            "is_pseudogene":
                self.is_pseudogene,

            "sequence_nt":
                self.sequence_nt,

            "region_length_nt":
                self.region_length_nt,

            "reference_coverage":
                self.reference_coverage,

            "is_small_region":
                self.is_small_region,

            "orfs":
                self.orfs or [],

            "num_orfs":
                self.num_orfs,

            "best_orf_length_nt":
                self.best_orf_length_nt,

            "mutation_evidence":
                self.mutation_evidence or {},
        })

        # Keep original TBLASTN protein alignments
        # for backwards compatibility.
        hsps = sorted(
            self.region_annotation.best_protein.hsps,
            key=lambda h: h.qstart
        )

        result["qseq"] = "".join(
            h.qseq for h in hsps
        )

        result["sseq"] = "".join(
            h.sseq for h in hsps
        )

        # These concatenated sequences are not
        # reconstructed full-length proteins.
        result["hsp_sequences_are_concatenated"] = True

        return result

    # --------------------------------------------------------
    # TSV SERIALIZATION
    # --------------------------------------------------------

    def to_tsv_dict(self) -> Dict[str, Any]:
        """
        Flat representation for final TSV output.
        """

        result = self.region_annotation.to_dict()

        result.update(self.disablements.to_dict())

        result.update({
            "classification_status":
                self.classification_status,

            "is_small_orf":
                self.is_small_orf,

            "is_pseudogene":
                self.is_pseudogene,

            "sequence_nt":
                self.sequence_nt,

            "region_length_nt":
                self.region_length_nt,

            "reference_coverage":
                self.reference_coverage,

            "is_small_region":
                self.is_small_region,

            "num_orfs":
                self.num_orfs,

            "best_orf_length_nt":
                self.best_orf_length_nt,

            # Store structured observations as JSON
            # inside a TSV cell.
            "mutation_evidence":
                __import__("json").dumps(
                    self.mutation_evidence or {}
                ),
        })

        hsps = sorted(
            self.region_annotation.best_protein.hsps,
            key=lambda h: h.qstart
        )

        result["qseq"] = "".join(
            h.qseq for h in hsps
        )

        result["sseq"] = "".join(
            h.sseq for h in hsps
        )

        return result

    # --------------------------------------------------------
    # DESERIALIZATION
    # --------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any]
    ) -> "PseudogeneAnnotation":
        """
        Reconstruct a candidate annotation from JSON.

        Preserves:
            - SSEARCH parent selection
            - TBLASTN HSPs
            - Genomic nucleotide sequence
            - ORF measurements
            - Mutation observations
            - Classification
        """

        return cls(
            region_annotation=RegionAnnotation.from_dict(d),

            disablements=DisablementCounts.from_dict(d),

            classification_status=d.get(
                "classification_status",
                "candidate"
            ),

            is_small_orf=d.get(
                "is_small_orf", False
            ),

            is_pseudogene=d.get(
                "is_pseudogene"
            ),

            sequence_nt=d.get(
                "sequence_nt"
            ),

            region_length_nt=d.get(
                "region_length_nt"
            ),

            reference_coverage=d.get(
                "reference_coverage"
            ),

            is_small_region=d.get(
                "is_small_region"
            ),

            orfs=d.get("orfs") or [],

            best_orf_length_nt=d.get(
                "best_orf_length_nt"
            ),

            mutation_evidence=d.get(
                "mutation_evidence"
            ) or {},
        )
