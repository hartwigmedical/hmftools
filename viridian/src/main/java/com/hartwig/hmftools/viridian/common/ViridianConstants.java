package com.hartwig.hmftools.viridian.common;

import java.util.List;
import java.util.Set;

public class ViridianConstants
{
    public static final String APP_NAME = "Viridian";

    // Output file suffixes.
    public static final String CANDIDATE_READ_FASTA_SUFFIX = ".viridian.reads.fasta";
    public static final String ALIGNED_READ_BAM_SUFFIX = ".viridian.reads.bam";
    public static final String CONTIG_INFO_TSV_SUFFIX = ".viridian.contig_info.tsv";
    public static final String PAIRWISE_MARGINS_TSV_SUFFIX = ".viridian.pairwise_margins.tsv";
    public static final String INTEGRATIONS_TSV_SUFFIX = ".viridian.integrations.tsv";

    // Genome partition size for multi-threaded candidate extraction.
    public static final int VIRAL_READ_EXTRACTION_PARTITION_SIZE = 1_000_000;

    // Contigs in the reference genome whose mapped reads are viral candidates.
    // Currently only applicable to v38 genome.
    public static final Set<String> VIRAL_REF_CONTIGS = Set.of("chrEBV");

    // Minimum soft-clip length for a host-mapped read to count as a candidate.
    public static final int VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT = 30;

    // Minimum BWA-MEM alignment score (-T), our only score floor. Set explicitly for visibility; matches the BWA default.
    public static final int VIRAL_READ_MIN_ALIGNMENT_SCORE_DEFAULT = 30;

    // Candidate reads submitted to BWA per alignment call (bounds memory usage).
    public static final int VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT = 100_000;

    // How softly a contested read's vote splits across strains: chance a base is right. Pessimistically low, as per-base
    // qualities are not retained.
    public static final double READ_VOTE_CORRECT_BASE_PROBABILITY = 0.25;

    // A clip counts as hanging over a contig end (a circular-genome artifact) when its bases would extend past the
    // boundary by more than this tolerance. Such alignments are dropped from all stats.
    public static final int VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE = 4;

    // An oncology group is present only if one of its contigs reaches this coverage fraction.
    public static final double VIRAL_CONTIG_COVERAGE_MIN = 0.1;

    // If any contig in an oncology group passes VIRAL_CONTIG_COVERAGE_MIN, then use this as the coverage threshold instead of that.
    // Ensures that a similar sibling contig which straddles the coverage threshold isn't harshly lost.
    public static final double VIRAL_CONTIG_COVERAGE_MIN_LOWER = 0.09;

    // Minimum read votes per base that is required to be covered.
    // Handles cases where the contig has some alignments but those reads strongly prefer a better matching contig
    // (since coverage is calculated based on ANY alignment, not necessarily a good one).
    public static final double VIRAL_CONTIG_VOTES_PER_BASE_MIN = 1.0;

    // Effectively a tolerance on the read vote share to ensure that a close runner up is also considered.
    // Otherwise you can have the top contig win by only an arbitrarily small number of read votes.
    public static final double REPRESENTATIVE_COMPARABLE_VOTE_RATIO = 0.9;

    // How many bases better does an alignment need to be to classify as diverging?
    // Lower is more sensitive to rival virus strains.
    public static final int REPRESENTATIVE_CHALLENGE_MARGIN_MIN = 5;

    // What fraction of the reads aligning to an oncology group must diverge to identify a possible rival virus strain?
    // The denominator is the distinct reads aligning anywhere in that oncology group.
    // Lower is more sensitive to rival virus strains.
    public static final double REPRESENTATIVE_CHALLENGE_READS_MIN = 0.1;

    // Winning margins reported per contig pair in the verbose output (for debugging and tuning only).
    public static final List<Integer> REPORTED_MARGINS = List.of(1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50);

    // Thresholds for insert sequence length to identify variants which could be viral integrations.
    // SGL length is reduced because usually SGL extension assembly is shorter.
    public static final int INTEGRATION_SGL_INSERT_LENGTH_MIN = 20;
    public static final int INTEGRATION_VARIANT_INSERT_LENGTH_MIN = 50;

    public static final int INTEGRATION_ALIGN_SCORE_MIN = 20;

    static
    {
        if(!(VIRAL_READ_EXTRACTION_PARTITION_SIZE > 0 && VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT > 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT > 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_READ_MIN_ALIGNMENT_SCORE_DEFAULT >= 0))
        {
            throw new IllegalStateException();
        }
        if(!(READ_VOTE_CORRECT_BASE_PROBABILITY > 0.0 && READ_VOTE_CORRECT_BASE_PROBABILITY < 1.0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE >= 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_CONTIG_COVERAGE_MIN > 0.0 && VIRAL_CONTIG_COVERAGE_MIN <= 1.0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_CONTIG_COVERAGE_MIN_LOWER > 0.0 && VIRAL_CONTIG_COVERAGE_MIN_LOWER <= VIRAL_CONTIG_COVERAGE_MIN))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_CONTIG_VOTES_PER_BASE_MIN > 0.0))
        {
            throw new IllegalStateException();
        }
        if(!(REPRESENTATIVE_COMPARABLE_VOTE_RATIO > 0.0 && REPRESENTATIVE_COMPARABLE_VOTE_RATIO <= 1.0))
        {
            throw new IllegalStateException();
        }
        if(!(REPRESENTATIVE_CHALLENGE_MARGIN_MIN >= 1))
        {
            throw new IllegalStateException();
        }
        if(!(REPRESENTATIVE_CHALLENGE_READS_MIN > 0.0 && REPRESENTATIVE_CHALLENGE_READS_MIN <= 1.0))
        {
            throw new IllegalStateException();
        }
        if(REPORTED_MARGINS.isEmpty())
        {
            throw new IllegalStateException();
        }
        if(!(INTEGRATION_SGL_INSERT_LENGTH_MIN > 0 && INTEGRATION_SGL_INSERT_LENGTH_MIN <= INTEGRATION_VARIANT_INSERT_LENGTH_MIN))
        {
            throw new IllegalStateException();
        }
        if(!(INTEGRATION_ALIGN_SCORE_MIN >= 0))
        {
            throw new IllegalStateException();
        }
    }
}
