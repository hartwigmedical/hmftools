package com.hartwig.hmftools.viridian.app;

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
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public record ViridianConfig(
        String sampleId,
        String tumorBam,
        String refGenomeFile,
        String viralRefFile,
        String viralRefInfoFile,
        String viralBwaIndexImage,
        @Nullable String bwaLibPath,
        @Nullable String esveeUnfilteredVcf,
        int alignmentBatchSize,
        int threads,
        String outputDir,
        @Nullable String outputId,
        boolean reuseReadsFasta,
        boolean reuseReadsBam,
        boolean verboseOutput
)
{
    private static final String CFG_VIRAL_REF_FILE = "viral_ref";
    private static final String DESC_VIRAL_REF_FILE = "Curated viral reference FASTA file";
    private static final String CFG_VIRAL_REF_INFO_FILE = "viral_ref_info";
    private static final String DESC_VIRAL_REF_INFO_FILE = "Viral reference info TSV (contig -> virus name + oncology group)";
    private static final String CFG_VIRAL_BWA_INDEX_IMAGE = "viral_bwa_index_image";
    private static final String DESC_VIRAL_BWA_INDEX_IMAGE = "Viral reference BWA-MEM index GATK image file";
    private static final String CFG_ESVEE_UNFILTERED_VCF = "esvee_unfiltered_vcf";
    private static final String DESC_ESVEE_UNFILTERED_VCF = "ESVEE caller unfiltered VCF";
    private static final String CFG_ALIGNMENT_BATCH_SIZE = "align_batch_size";
    private static final String DESC_ALIGNMENT_BATCH_SIZE = "Candidate reads submitted to BWA per alignment call";
    private static final String CFG_REUSE_READS_FASTA = "reuse_reads_fasta";
    private static final String DESC_REUSE_READS_FASTA =
            "Dev: skip read extraction and reuse the candidate read FASTA already at the output location";
    private static final String CFG_REUSE_READS_BAM = "reuse_reads_bam";
    private static final String DESC_REUSE_READS_BAM =
            "Dev: skip read extraction and alignment, and reuse the aligned read BAM already at the output location";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information which may be useful for debugging";

    public static ViridianConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String viralRefFile = configBuilder.getValue(CFG_VIRAL_REF_FILE);
        return new ViridianConfig(
                configBuilder.getValue(SAMPLE),
                configBuilder.getValue(TUMOR_BAM),
                configBuilder.getValue(REF_GENOME),
                viralRefFile,
                configBuilder.getValue(CFG_VIRAL_REF_INFO_FILE),
                configBuilder.getValue(CFG_VIRAL_BWA_INDEX_IMAGE, viralRefFile + ".img"),
                configBuilder.getValue(BWA_LIB_PATH),
                configBuilder.getValue(CFG_ESVEE_UNFILTERED_VCF),
                configBuilder.getInteger(CFG_ALIGNMENT_BATCH_SIZE),
                parseThreads(configBuilder),
                parseOutputDir(configBuilder),
                configBuilder.getValue(OUTPUT_ID),
                configBuilder.hasFlag(CFG_REUSE_READS_FASTA),
                configBuilder.hasFlag(CFG_REUSE_READS_BAM),
                configBuilder.hasFlag(CFG_VERBOSE_OUTPUT)
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
        configBuilder.addPath(CFG_ESVEE_UNFILTERED_VCF, false, DESC_ESVEE_UNFILTERED_VCF);

        configBuilder.addInteger(CFG_ALIGNMENT_BATCH_SIZE, DESC_ALIGNMENT_BATCH_SIZE, VIRAL_READ_ALIGNMENT_BATCH_SIZE_DEFAULT);
        configBuilder.addFlag(CFG_REUSE_READS_FASTA, DESC_REUSE_READS_FASTA);
        configBuilder.addFlag(CFG_REUSE_READS_BAM, DESC_REUSE_READS_BAM);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        addThreadOptions(configBuilder);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        addOutputId(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
