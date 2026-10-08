package com.hartwig.hmftools.viridian.app;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_TUMOR_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BAM_METRICS_TUMOR_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputId;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT;

import java.util.Arrays;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.detection.extract.CandidateReadSource;

import org.jetbrains.annotations.Nullable;

public record ViridianConfig(
        String sampleId,
        String tumorBam,
        String refGenomeFile,
        String virusRefFile,
        String virusRefInfoFile,
        String virusBwaIndexImage,
        String oncologyGroupInfoFile,
        @Nullable String bwaLibPath,
        String esveeUnfilteredVcf,
        String purpleDir,
        String bamMetricsTumorDir,
        int alignmentBatchSize,
        int threads,
        String outputDir,
        @Nullable String outputId,
        CandidateReadSource candidateReadSource,
        boolean kmerFilterEnabled,
        boolean verboseOutput
)
{
    private static final String CFG_VIRUS_REF_FILE = "virus_ref";
    private static final String DESC_VIRUS_REF_FILE = "Curated virus reference FASTA file";
    private static final String CFG_VIRUS_REF_INFO_FILE = "virus_ref_info";
    private static final String DESC_VIRUS_REF_INFO_FILE = "Virus reference info TSV (contig -> virus name + oncology group)";
    private static final String CFG_VIRUS_BWA_INDEX_IMAGE = "virus_bwa_index_image";
    private static final String DESC_VIRUS_BWA_INDEX_IMAGE = "Virus reference BWA-MEM index GATK image file";
    private static final String CFG_ESVEE_UNFILTERED_VCF = "esvee_unfiltered_vcf";
    private static final String DESC_ESVEE_UNFILTERED_VCF = "ESVEE caller unfiltered VCF";
    private static final String CFG_ONCOLOGY_GROUP_INFO_FILE = "oncology_group_info";
    private static final String DESC_ONCOLOGY_GROUP_INFO_FILE = "Virus oncology group info TSV";
    private static final String CFG_ALIGNMENT_BATCH_SIZE = "align_batch_size";
    private static final String DESC_ALIGNMENT_BATCH_SIZE = "Candidate reads submitted to BWA per alignment call";
    static final String CFG_CANDIDATE_READ_SOURCE = "candidate_read_source";
    private static final String CANDIDATE_READ_SOURCE_VALUES =
            Arrays.stream(CandidateReadSource.values()).map(Enum::name).collect(Collectors.joining(", "));
    private static final String DESC_CANDIDATE_READ_SOURCE = "Candidate read source: " + CANDIDATE_READ_SOURCE_VALUES;
    static final String CFG_KMER_FILTER = "kmer_filter";
    private static final String DESC_KMER_FILTER = "Filter candidate viral reads with a kmer index";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information which may be useful for debugging";

    public static ViridianConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String virusRefFile = configBuilder.getValue(CFG_VIRUS_REF_FILE);
        return new ViridianConfig(
                configBuilder.getValue(SAMPLE),
                configBuilder.getValue(TUMOR_BAM),
                configBuilder.getValue(REF_GENOME),
                virusRefFile,
                configBuilder.getValue(CFG_VIRUS_REF_INFO_FILE),
                configBuilder.getValue(CFG_VIRUS_BWA_INDEX_IMAGE, virusRefFile + ".img"),
                configBuilder.getValue(CFG_ONCOLOGY_GROUP_INFO_FILE),
                configBuilder.getValue(BWA_LIB_PATH),
                configBuilder.getValue(CFG_ESVEE_UNFILTERED_VCF),
                configBuilder.getValue(PURPLE_DIR_CFG),
                configBuilder.getValue(BAM_METRICS_TUMOR_DIR_CFG),
                configBuilder.getInteger(CFG_ALIGNMENT_BATCH_SIZE),
                parseThreads(configBuilder),
                parseOutputDir(configBuilder),
                configBuilder.getValue(OUTPUT_ID),
                parseCandidateReadSource(configBuilder.getValue(CFG_CANDIDATE_READ_SOURCE)),
                configBuilder.hasFlag(CFG_KMER_FILTER),
                configBuilder.hasFlag(CFG_VERBOSE_OUTPUT)
        );
    }

    private static CandidateReadSource parseCandidateReadSource(String value)
    {
        try
        {
            return CandidateReadSource.valueOf(value);
        }
        catch(IllegalArgumentException e)
        {
            throw new UserInputError(format(
                    "Invalid %s: %s, expected one of: %s", CFG_CANDIDATE_READ_SOURCE, value, CANDIDATE_READ_SOURCE_VALUES));
        }
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(TUMOR_BAM, true, TUMOR_BAM_DESC);
        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(CFG_VIRUS_REF_FILE, true, DESC_VIRUS_REF_FILE);
        configBuilder.addPath(CFG_VIRUS_REF_INFO_FILE, true, DESC_VIRUS_REF_INFO_FILE);
        configBuilder.addPath(CFG_VIRUS_BWA_INDEX_IMAGE, false, DESC_VIRUS_BWA_INDEX_IMAGE);
        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);
        configBuilder.addPath(CFG_ESVEE_UNFILTERED_VCF, true, DESC_ESVEE_UNFILTERED_VCF);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(BAM_METRICS_TUMOR_DIR_CFG, true, BAM_METRICS_TUMOR_DIR_DESC);
        configBuilder.addPath(CFG_ONCOLOGY_GROUP_INFO_FILE, true, DESC_ONCOLOGY_GROUP_INFO_FILE);

        configBuilder.addInteger(CFG_ALIGNMENT_BATCH_SIZE, DESC_ALIGNMENT_BATCH_SIZE, VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT);
        configBuilder.addConfigItem(
                CFG_CANDIDATE_READ_SOURCE, false, DESC_CANDIDATE_READ_SOURCE, CandidateReadSource.ALL.toString());
        configBuilder.addFlag(CFG_KMER_FILTER, DESC_KMER_FILTER);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        addThreadOptions(configBuilder);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        addOutputId(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
