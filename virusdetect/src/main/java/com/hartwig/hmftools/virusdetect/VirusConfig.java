package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
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
import static com.hartwig.hmftools.virusdetect.VirusConstants.ALIGNMENT_BATCH_SIZE_DEFAULT;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public record VirusConfig(
        String sampleId,
        String tumorBam,
        String refGenomeFile,
        String viralRefFile,
        String viralRefInfoFile,
        String viralBwaIndexImage,
        @Nullable String bwaLibPath,
        int alignmentBatchSize,
        int threads,
        String outputDir,
        @Nullable String outputId
)
{
    private static final String CFG_VIRAL_REF_FILE = "viral_ref";
    private static final String DESC_VIRAL_REF_FILE = "Curated viral reference FASTA file";
    private static final String CFG_VIRAL_REF_INFO_FILE = "viral_ref_info";
    private static final String DESC_VIRAL_REF_INFO_FILE = "Viral reference info TSV (contig -> virus name + oncology group)";
    private static final String CFG_VIRAL_BWA_INDEX_IMAGE = "viral_bwa_index_image";
    private static final String DESC_VIRAL_BWA_INDEX_IMAGE = "Viral reference BWA-MEM index GATK image file";
    private static final String CFG_ALIGNMENT_BATCH_SIZE = "align_batch_size";
    private static final String DESC_ALIGNMENT_BATCH_SIZE = "Candidate reads submitted to BWA per alignment call";

    public static VirusConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String viralRefFile = configBuilder.getValue(CFG_VIRAL_REF_FILE);
        return new VirusConfig(
                configBuilder.getValue(SAMPLE),
                configBuilder.getValue(TUMOR_BAM),
                configBuilder.getValue(REF_GENOME),
                viralRefFile,
                configBuilder.getValue(CFG_VIRAL_REF_INFO_FILE),
                configBuilder.getValue(CFG_VIRAL_BWA_INDEX_IMAGE, viralRefFile + ".img"),
                configBuilder.getValue(BWA_LIB_PATH),
                configBuilder.getInteger(CFG_ALIGNMENT_BATCH_SIZE),
                parseThreads(configBuilder),
                parseOutputDir(configBuilder),
                configBuilder.getValue(OUTPUT_ID)
        );
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(TUMOR_BAM, true, TUMOR_BAM_DESC);
        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(CFG_VIRAL_REF_FILE, true, DESC_VIRAL_REF_FILE);
        configBuilder.addPath(CFG_VIRAL_REF_INFO_FILE, true, DESC_VIRAL_REF_INFO_FILE);
        configBuilder.addPath(CFG_VIRAL_BWA_INDEX_IMAGE, false, DESC_VIRAL_BWA_INDEX_IMAGE);
        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);

        configBuilder.addInteger(CFG_ALIGNMENT_BATCH_SIZE, DESC_ALIGNMENT_BATCH_SIZE, ALIGNMENT_BATCH_SIZE_DEFAULT);

        addThreadOptions(configBuilder);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        addOutputId(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
