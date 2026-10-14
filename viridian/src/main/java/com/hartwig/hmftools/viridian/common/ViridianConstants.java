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
    public static final String GENOME_SUPPORT_TSV_SUFFIX = ".genome_support.tsv";
    public static final String PAIRWISE_MARGINS_TSV_SUFFIX = ".pairwise_margins.tsv";
    public static final String VIRUS_DETECTION_TSV_SUFFIX = ".virus_detection.tsv";
    public static final String INTEGRATION_CANDIDATES_TSV_SUFFIX = ".integration_candidates.tsv";
    public static final String VIRUS_REPORT_TSV_SUFFIX = ".virus_report.tsv";
    public static final String RUN_MANIFEST_TSV_SUFFIX = ".manifest.tsv";

    // Genome partition size for multi-threaded candidate extraction.
    // This value was chosen approximately based on observed highest performance.
    public static final int VIRAL_READ_EXTRACTION_PARTITION_SIZE = 1_000_000;

    // Contigs in the reference genome whose mapped reads are viral candidates.
    public static final Set<String> DEFAULT_VIRUS_DECOY_CONTIGS = Set.of(
            // Epstein-Barr virus, in GRCh38 analysis sets.
            "chrEBV",
            // Epstein-Barr virus, in GRCh37 hs37d5 and GATK b37.
            "NC_007605",
            // The remaining contigs are from GDC GRCh38.d1.vd1.
            // Human cytomegalovirus.
            "CMV",
            // Hepatitis B virus.
            "HBV",
            // Hepatitis C virus.
            "HCV-1", "HCV-2",
            // Human immunodeficiency virus.
            "HIV-1", "HIV-2",
            // Kaposi sarcoma-associated herpesvirus.
            "KSHV",
            // Human T-lymphotropic virus 1.
            "HTLV-1",
            // Merkel cell polyomavirus.
            "MCV",
            // Simian virus 40.
            "SV40",
            // Human papillomavirus types.
            "HPV1", "HPV2", "HPV3", "HPV4", "HPV5", "HPV6", "HPV7", "HPV8", "HPV9", "HPV10", "HPV11", "HPV12", "HPV13", "HPV14",
            "HPV15", "HPV16", "HPV17", "HPV18", "HPV19", "HPV20", "HPV21", "HPV22", "HPV23", "HPV24", "HPV25", "HPV26", "HPV27",
            "HPV28", "HPV29", "HPV30", "HPV31", "HPV32", "HPV33", "HPV34", "HPV35", "HPV36", "HPV37", "HPV38", "HPV39", "HPV40",
            "HPV41", "HPV42", "HPV43", "HPV44", "HPV45", "HPV47", "HPV48", "HPV49", "HPV50", "HPV51", "HPV52", "HPV53", "HPV54",
            "HPV56", "HPV57", "HPV58", "HPV59", "HPV60", "HPV61", "HPV62", "HPV63", "HPV65", "HPV66", "HPV67", "HPV68", "HPV69",
            "HPV70", "HPV71", "HPV72", "HPV73", "HPV74", "HPV75", "HPV76", "HPV77", "HPV78", "HPV80", "HPV81", "HPV82", "HPV83",
            "HPV84", "HPV85", "HPV86", "HPV87", "HPV88", "HPV89", "HPV90", "HPV91", "HPV92", "HPV93", "HPV94", "HPV95", "HPV96",
            "HPV97", "HPV98", "HPV99", "HPV100", "HPV101", "HPV102", "HPV103", "HPV104", "HPV105", "HPV106", "HPV107", "HPV108",
            "HPV109", "HPV110", "HPV111", "HPV112", "HPV113", "HPV114", "HPV115", "HPV116", "HPV117", "HPV118", "HPV119",
            "HPV120", "HPV121", "HPV122", "HPV123", "HPV124", "HPV125", "HPV126", "HPV127", "HPV128", "HPV129", "HPV130",
            "HPV131", "HPV132", "HPV133", "HPV134", "HPV135", "HPV136", "HPV137", "HPV138", "HPV139", "HPV140", "HPV141",
            "HPV142", "HPV143", "HPV144", "HPV145", "HPV146", "HPV147", "HPV148", "HPV149", "HPV150", "HPV151", "HPV152",
            "HPV153", "HPV154", "HPV155", "HPV156", "HPV159", "HPV160", "HPV161", "HPV162", "HPV163", "HPV164", "HPV165",
            "HPV166", "HPV167", "HPV168", "HPV169", "HPV170", "HPV171", "HPV172", "HPV173", "HPV174", "HPV175", "HPV178",
            "HPV179", "HPV180", "HPV184", "HPV197", "HPV199",
            "HPV-mCG2", "HPV-mCG3", "HPV-mCH2", "HPV-mFD1", "HPV-mFD2", "HPV-mFS1", "HPV-mFi864", "HPV-mKC5", "HPV-mKN1",
            "HPV-mKN2", "HPV-mKN3", "HPV-mL55", "HPV-mRTRX7", "HPV-mSD2");

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

    // Affects attribution of read votes between virus genomes.
    // This value is pessimistic - assumes the correct base is completely unknown.
    public static final double READ_VOTE_CORRECT_BASE_PROBABILITY = 0.25;

    // A clip counts as hanging over a genome end (a circular-genome artifact) when its bases would extend past the
    // boundary by more than this tolerance. Such alignments are dropped from all further processing.
    public static final int VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE = 4;

    // A genome is a candidate only if some genome in its oncology group reaches this coverage fraction.
    // This value is very conservative. Hard to argue a lower coverage wouldn't be spurious.
    public static final double VIRUS_GENOME_COVERAGE_MIN = 0.1;

    // If any genome in an oncology group passes VIRUS_GENOME_COVERAGE_MIN, then use this as the coverage threshold instead of that.
    // Ensures that a similar sibling genome which straddles the coverage threshold isn't harshly lost.
    public static final double VIRUS_GENOME_COVERAGE_MIN_LOWER = VIRUS_GENOME_COVERAGE_MIN * 0.9;

    // Minimum read votes per base that is required to be covered.
    // Handles cases where the genome has some alignments but those reads strongly prefer a better matching genome
    // (since coverage is calculated based on ANY alignment, not necessarily a good one).
    public static final double VIRUS_GENOME_VOTES_PER_BASE_MIN = 1.0;

    // Effectively a tolerance on the read vote share to ensure that a close runner up is also considered.
    // Otherwise, you can have the top genome win by only an arbitrarily small number of read votes.
    public static final double REPRESENTATIVE_COMPARABLE_VOTE_RATIO = 0.9;

    // How many bases better does an alignment need to be to classify as diverging?
    // Lower is more sensitive to rival virus strains.
    public static final int REPRESENTATIVE_CHALLENGE_MARGIN_MIN = 5;

    // What fraction of the reads aligning to an oncology group must diverge to identify a possible rival virus strain?
    // The denominator is the distinct reads aligning anywhere in that oncology group.
    // Lower is more sensitive to rival virus strains.
    public static final double REPRESENTATIVE_CHALLENGE_READS_MIN = 0.1;

    // Winning margins reported per genome pair in the verbose output (for debugging and tuning only).
    public static final List<Integer> REPORTED_MARGINS = List.of(1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50);

    // Minimum insert sequence length for an SV to be a viral integration candidate.
    // This is a conservative value. ESVEE doesn't call events shorter than 32b anyway.
    public static final int INTEGRATION_INSERT_LENGTH_MIN = 30;

    // Minimum alignment score for host insert sequences aligned to virus genomes.
    // This is a conservative value, validated with observed concordance with VirusBreakend.
    public static final int INTEGRATION_ALIGN_SCORE_MIN = 30;

    // Minimum aligner score per aligned base. Used to filter out alignments which are long but poor similarity.
    public static final double INTEGRATION_VIRAL_ALIGN_SCORE_PER_BASE_MIN = 0.7;

    // Minimum integrations for a virus to be reported.
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
        if(!(VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE >= 0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRUS_GENOME_COVERAGE_MIN > 0.0 && VIRUS_GENOME_COVERAGE_MIN <= 1.0))
        {
            throw new IllegalStateException();
        }
        if(!(VIRUS_GENOME_COVERAGE_MIN_LOWER > 0.0 && VIRUS_GENOME_COVERAGE_MIN_LOWER <= VIRUS_GENOME_COVERAGE_MIN))
        {
            throw new IllegalStateException();
        }
        if(!(VIRUS_GENOME_VOTES_PER_BASE_MIN > 0.0))
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
