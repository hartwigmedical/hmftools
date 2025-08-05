package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.mappability.ProbeQualityProfile.CFG_PROBE_QUALITY_FILE;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariantsConfig;

import org.jetbrains.annotations.Nullable;

public record PanelBuilderConfig(
        String refGenomeFile,
        String ensemblDir,
        String probeQualityProfileFile,
        @Nullable String amberSitesFile,
        String bwaIndexImageFile,
        @Nullable String bwaLibPath,
        @Nullable String genesFile,
        boolean includeCdr3,
        @Nullable String customRegionsFile,
        @Nullable SampleVariantsConfig sampleVariants,
        @Nullable String outputId,
        String outputDir,
        boolean verboseOutput,
        int threads
)
{
    private static final String CFG_BWA_INDEX_IMAGE_FILE = "bwa_index_image";
    private static final String DESC_BWA_INDEX_IMAGE_FILE = "Reference genome BWA-MEM index GATK image file";
    private static final String CFG_AMBER_SITES_FILE = "amber_sites";
    private static final String DESC_AMBER_SITES_FILE = "Amber het sites file";
    private static final String CFG_TARGET_GENES_FILE = "target_genes";
    private static final String DESC_TARGET_GENES_FILE = "Gene and transcript name file";
    private static final String CFG_INCLUDE_CDR3 = "cdr3";
    private static final String DESC_INCLUDE_CDR3 = "Include fixed CDR3 panel probes";
    private static final String CFG_CUSTOM_REGIONS_FILE = "custom_regions";
    private static final String DESC_CUSTOM_REGIONS_FILE = "Custom region file";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information useful for debugging";

    public static PanelBuilderConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String refGenomePath = configBuilder.getValue(REF_GENOME);
        return new PanelBuilderConfig(
                refGenomePath,
                configBuilder.getValue(ENSEMBL_DATA_DIR),
                configBuilder.getValue(CFG_PROBE_QUALITY_FILE),
                configBuilder.getValue(CFG_AMBER_SITES_FILE),
                configBuilder.getValue(CFG_BWA_INDEX_IMAGE_FILE, refGenomePath + ".img"),
                configBuilder.getValue(BWA_LIB_PATH),
                configBuilder.getValue(CFG_TARGET_GENES_FILE),
                configBuilder.hasFlag(CFG_INCLUDE_CDR3),
                configBuilder.getValue(CFG_CUSTOM_REGIONS_FILE),
                SampleVariantsConfig.fromConfigBuilder(configBuilder),
                configBuilder.getValue(OUTPUT_ID),
                parseOutputDir(configBuilder),
                configBuilder.hasFlag(CFG_VERBOSE_OUTPUT),
                parseThreads(configBuilder)
        );
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        addEnsemblDir(configBuilder, true);
        ProbeQualityProfile.registerConfig(configBuilder);
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, DESC_AMBER_SITES_FILE);

        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);
        configBuilder.addPath(CFG_BWA_INDEX_IMAGE_FILE, false, DESC_BWA_INDEX_IMAGE_FILE);

        configBuilder.addPath(CFG_TARGET_GENES_FILE, false, DESC_TARGET_GENES_FILE);
        configBuilder.addFlag(CFG_INCLUDE_CDR3, DESC_INCLUDE_CDR3);
        configBuilder.addPath(CFG_CUSTOM_REGIONS_FILE, false, DESC_CUSTOM_REGIONS_FILE);

        addOutputOptions(configBuilder);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        addThreadOptions(configBuilder);

        addLoggingOptions(configBuilder);

        SampleVariantsConfig.registerConfig(configBuilder);
    }
}
