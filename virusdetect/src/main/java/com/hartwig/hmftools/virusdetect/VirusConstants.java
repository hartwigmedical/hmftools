package com.hartwig.hmftools.virusdetect;

import java.util.List;
import java.util.Set;

public class VirusConstants
{
    public static final String APP_NAME = "VirusDetect";

    // Output file suffixes.
    public static final String CANDIDATE_FASTA_SUFFIX = ".virus.candidates.fasta";
    public static final String ALIGNED_BAM_SUFFIX = ".virus.aligned.bam";
    public static final String CONTIG_INFO_TSV_SUFFIX = ".virus.contig_info.tsv";
    public static final String PAIRWISE_MARGINS_TSV_SUFFIX = ".virus.pairwise_margins.tsv";

    // Genome partition size for multi-threaded candidate extraction.
    public static final int EXTRACTION_PARTITION_SIZE = 1_000_000;

    // Contigs in the reference genome whose mapped reads are viral candidates (v38 naming).
    public static final Set<String> DECOY_CONTIGS = Set.of("chrEBV");

    // Minimum soft-clip length for a host-mapped read to count as a candidate.
    public static final int MIN_SOFT_CLIP_BASES_DEFAULT = 30;

    // Minimum BWA-MEM alignment score (-T), our only score floor. Set explicitly for visibility; matches the BWA default.
    public static final int MIN_ALIGNMENT_SCORE_DEFAULT = 30;

    // Candidate reads submitted to BWA per alignment call (bounds memory usage).
    public static final int ALIGNMENT_BATCH_SIZE_DEFAULT = 100_000;

    // How softly a contested read's vote splits across strains: chance a base is right. Pessimistically low, as per-base
    // qualities are not retained.
    public static final double VOTE_CORRECT_BASE_PROBABILITY = 0.25;

    // A clip counts as hanging over a contig end (a circular-genome artifact) when its bases would extend past the
    // boundary by more than this tolerance. Such alignments are dropped from all stats.
    public static final int ORIGIN_CLIP_TOLERANCE = 4;

    // An oncology group is present only if one of its contigs reaches this coverage fraction.
    public static final double MIN_COVERAGE = 0.1;

    // If any contig in an oncology group passes MIN_COVERAGE, then use this as the coverage threshold instead of MIN_COVERAGE.
    // Ensures that a similar sibling contig which straddles the coverage threshold isn't harshly lost.
    public static final double MIN_COVERAGE_LOWER = 0.09;

    // Minimum read votes per base that is required to be covered.
    // Handles cases where the contig has some alignments but those reads strongly prefer a better matching contig
    // (since coverage is calculated based on ANY alignment, not necessarily a good one).
    public static final double MIN_VOTES_PER_BASE = 1.0;

    // Effectively a tolerance on the read vote share to ensure that a close runner up is also considered.
    // Otherwise you can have the top contig win by only an arbitrarily small number of read votes.
    public static final double COMPARABLE_VOTE_RATIO = 0.9;

    // How many bases better does an alignment need to be to classify as diverging?
    // Lower is more sensitive to rival virus strains.
    public static final int MIN_CHALLENGE_MARGIN = 5;

    // What fraction of the reads aligning to an oncology group must diverge to identify a possible rival virus strain?
    // The denominator is the distinct reads aligning anywhere in that oncology group.
    // Lower is more sensitive to rival virus strains.
    public static final double MIN_CHALLENGE_READS = 0.1;

    // Winning margins reported per contig pair in the verbose output (for debugging and tuning only).
    public static final List<Integer> REPORTED_MARGINS = List.of(1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50);
}
