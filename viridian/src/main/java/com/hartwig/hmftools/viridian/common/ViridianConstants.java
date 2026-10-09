package com.hartwig.hmftools.viridian.common;

import java.util.List;
import java.util.Set;

public class ViridianConstants
{
    public static final String APP_NAME = "Viridian";

    // Output file suffixes. The app name is inserted where these build the full file name.
    public static final String CANDIDATES_FASTA_SUFFIX = ".candidates.fasta";
    public static final String ALL_ALIGNMENTS_BAM_SUFFIX = ".all.bam";
    public static final String REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX = ".representative.bam";
    public static final String CONTIG_SUPPORT_TSV_SUFFIX = ".contig_support.tsv";
    public static final String PAIRWISE_MARGINS_TSV_SUFFIX = ".pairwise_margins.tsv";
    public static final String VIRUS_DETECTION_TSV_SUFFIX = ".virus_detection.tsv";
    public static final String INTEGRATIONS_TSV_SUFFIX = ".integrations.tsv";
    public static final String VIRUS_REPORT_TSV_SUFFIX = ".virus_report.tsv";
    public static final String RUN_MANIFEST_TSV_SUFFIX = ".manifest.tsv";

    // Genome partition size for multi-threaded candidate extraction.
    // This value was chosen approximately based on observed highest performance.
    public static final int VIRAL_READ_EXTRACTION_PARTITION_SIZE = 1_000_000;

    // Contigs in the reference genome whose mapped reads are viral candidates.
    // Currently only applicable to v38 genome.
    public static final Set<String> VIRUS_REF_CONTIGS = Set.of("chrEBV");

    // Minimum soft-clip length for a host-mapped read to count as a candidate.
    public static final int VIRAL_READ_MIN_SOFT_CLIP_BASES = 30;

    // Reads must share an exact match of at least this length to a virus genome to be considered.
    // Should match the BWA-MEM seed size.
    public static final int VIRAL_KMER_LENGTH = 19;

    public static final int VIRAL_KMER_BLOOM_BITS = 1 << 27;

    public static final int VIRAL_KMER_BLOOM_HASHES = 3;

    // Minimum alignment score for reads aligning to virus genomes.
    public static final int VIRAL_READ_MIN_ALIGN_SCORE = 30;

    // Candidate reads submitted to BWA per alignment call (bounds memory usage).
    public static final int VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT = 100_000;

    // Affects attribution of read votes between contigs.
    // This value is pessimistic - assumes the correct base is completely unknown.
    public static final double READ_VOTE_CORRECT_BASE_PROBABILITY = 0.25;

    // A clip counts as hanging over a contig end (a circular-genome artifact) when its bases would extend past the
    // boundary by more than this tolerance. Such alignments are dropped from all further processing.
    public static final int VIRAL_CONTIG_ORIGIN_CLIP_TOLERANCE = 4;

    // A contig is a candidate only if some contig in its oncology group reaches this coverage fraction.
    // This value is very conservative. Hard to argue a lower coverage wouldn't be spurious.
    public static final double VIRAL_CONTIG_COVERAGE_MIN = 0.1;

    // If any contig in an oncology group passes VIRAL_CONTIG_COVERAGE_MIN, then use this as the coverage threshold instead of that.
    // Ensures that a similar sibling contig which straddles the coverage threshold isn't harshly lost.
    public static final double VIRAL_CONTIG_COVERAGE_MIN_LOWER = VIRAL_CONTIG_COVERAGE_MIN * 0.9;

    // Minimum read votes per base that is required to be covered.
    // Handles cases where the contig has some alignments but those reads strongly prefer a better matching contig
    // (since coverage is calculated based on ANY alignment, not necessarily a good one).
    public static final double VIRAL_CONTIG_VOTES_PER_BASE_MIN = 1.0;

    // Effectively a tolerance on the read vote share to ensure that a close runner up is also considered.
    // Otherwise, you can have the top contig win by only an arbitrarily small number of read votes.
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

    // Minimum insert sequence length for an SV to be a viral integration candidate.
    // This is a conservative value. ESVEE doesn't call events shorter than 32b anyway.
    public static final int INTEGRATION_INSERT_LENGTH_MIN = 30;

    // Minimum alignment score for host insert sequences aligned to virus genomes.
    // This is a conservative value, validated with observed concordance with VirusBreakend.
    public static final int INTEGRATION_ALIGN_SCORE_MIN = 30;

    // Minimum aligner score per aligned base. Used to filter out alignments which are long but poor similarity.
    public static final double INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN = 0.7;

    // Minimum plausible integrations for a virus to be reported.
    // This value is a starting point based on VirusInterpreter reporting logic.
    public static final int REPORTED_INTEGRATIONS_MIN = 1;

    // Minimum virus genome copies per tumor cell for a virus to be reported.
    // This value is a starting point based on VirusInterpreter reporting logic.
    public static final double REPORTED_COPIES_PER_CELL_MIN = 0.5;

    static
    {
        if(!(VIRAL_READ_EXTRACTION_PARTITION_SIZE > 0 && VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT > 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_READ_MIN_SOFT_CLIP_BASES > 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_KMER_LENGTH >= 19 && VIRAL_KMER_LENGTH <= 31))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_KMER_BLOOM_BITS >= 64 && Integer.bitCount(VIRAL_KMER_BLOOM_BITS) == 1))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_KMER_BLOOM_HASHES >= 1))
        {
            throw new IllegalStateException();
        }
        if(!(VIRAL_READ_MIN_ALIGN_SCORE >= 19))
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
        if(!(INTEGRATION_ALIGN_SCORE_MIN >= 19))
        {
            throw new IllegalStateException();
        }
        if(!(INTEGRATION_INSERT_LENGTH_MIN > 0 && INTEGRATION_INSERT_LENGTH_MIN <= INTEGRATION_ALIGN_SCORE_MIN))
        {
            throw new IllegalStateException();
        }
        if(!(INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN > 0.0 && INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN <= 1.0))
        {
            throw new IllegalStateException();
        }
        if(!(REPORTED_INTEGRATIONS_MIN >= 1))
        {
            throw new IllegalStateException();
        }
        if(!(REPORTED_COPIES_PER_CELL_MIN > 0.0))
        {
            throw new IllegalStateException();
        }
    }
}
