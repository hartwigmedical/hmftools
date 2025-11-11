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
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputId;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_RESOLUTION_KB_DEFAULT;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariantsConfig;

import org.jetbrains.annotations.Nullable;

public record PanelBuilderConfig(
        String refGenomeFile,
        @Nullable String ensemblDir,
        String probeQualityProfileFile,
        String bwaIndexImageFile,
        @Nullable String bwaLibPath,
        @Nullable String genesFile,
        boolean includeCnBackbone,
        @Nullable String amberSitesFile,
        int cnBackboneResolution,
        boolean includeCdr3,
        @Nullable String customRegionsFile,
        @Nullable String customSvsFile,
        @Nullable SampleVariantsConfig sampleVariants,
        int threads,
        @Nullable String outputId,
        String outputDir,
        boolean verboseOutput
)
{
    private static final String CFG_BWA_INDEX_IMAGE_FILE = "bwa_index_image";
    private static final String DESC_BWA_INDEX_IMAGE_FILE = "Reference genome BWA-MEM index GATK image file";
    private static final String CFG_INCLUDE_CN_BACKBONE = "cn_backbone";
    private static final String DESC_INCLUDE_CN_BACKBONE = "Include copy number backbone probes";
    private static final String CFG_AMBER_SITES_FILE = "amber_sites";
    private static final String DESC_AMBER_SITES_FILE = "Amber heterozygous sites TSV file";
    private static final String CFG_CN_BACKBONE_RESOLUTION = "cn_backbone_res_kb";
    private static final String DESC_CN_BACKBONE_RESOLUTION = "Approximate spacing between copy number backbone probes, in kb";
    private static final String CFG_TARGET_GENES_FILE = "genes";
    private static final String DESC_TARGET_GENES_FILE = "Gene options and transcript TSV file";
    private static final String CFG_INCLUDE_CDR3 = "cdr3";
    private static final String DESC_INCLUDE_CDR3 = "Include fixed CDR3 panel probes";
    private static final String CFG_CUSTOM_REGIONS_FILE = "custom_regions";
    private static final String DESC_CUSTOM_REGIONS_FILE = "Custom regions TSV file";
    private static final String CFG_CUSTOM_SVS_FILE = "custom_svs";
    private static final String DESC_CUSTOM_SV_FILE = "Custom structural variants TSV file";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information which may be useful for debugging";

    public static PanelBuilderConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String refGenomePath = configBuilder.getValue(REF_GENOME);
        return new PanelBuilderConfig(
                refGenomePath,
                configBuilder.getValue(ENSEMBL_DATA_DIR),
                configBuilder.getValue(CFG_PROBE_QUALITY_FILE),
                configBuilder.getValue(CFG_BWA_INDEX_IMAGE_FILE, refGenomePath + ".img"),
                configBuilder.getValue(BWA_LIB_PATH),
                configBuilder.getValue(CFG_TARGET_GENES_FILE),
                configBuilder.hasFlag(CFG_INCLUDE_CN_BACKBONE),
                configBuilder.getValue(CFG_AMBER_SITES_FILE),
                configBuilder.getInteger(CFG_CN_BACKBONE_RESOLUTION) * 1000,
                configBuilder.hasFlag(CFG_INCLUDE_CDR3),
                configBuilder.getValue(CFG_CUSTOM_REGIONS_FILE),
                configBuilder.getValue(CFG_CUSTOM_SVS_FILE),
                SampleVariantsConfig.fromConfigBuilder(configBuilder),
                parseThreads(configBuilder),
                configBuilder.getValue(OUTPUT_ID),
                parseOutputDir(configBuilder),
                configBuilder.hasFlag(CFG_VERBOSE_OUTPUT)
        );
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        addEnsemblDir(configBuilder, false);
        ProbeQualityProfile.registerConfig(configBuilder);
        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);
        configBuilder.addPath(CFG_BWA_INDEX_IMAGE_FILE, false, DESC_BWA_INDEX_IMAGE_FILE);

        configBuilder.addFlag(CFG_INCLUDE_CN_BACKBONE, DESC_INCLUDE_CN_BACKBONE);
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, DESC_AMBER_SITES_FILE);
        configBuilder.addInteger(CFG_CN_BACKBONE_RESOLUTION, DESC_CN_BACKBONE_RESOLUTION, CN_BACKBONE_RESOLUTION_KB_DEFAULT);
        configBuilder.addPath(CFG_TARGET_GENES_FILE, false, DESC_TARGET_GENES_FILE);
        configBuilder.addFlag(CFG_INCLUDE_CDR3, DESC_INCLUDE_CDR3);
        configBuilder.addPath(CFG_CUSTOM_REGIONS_FILE, false, DESC_CUSTOM_REGIONS_FILE);
        configBuilder.addPath(CFG_CUSTOM_SVS_FILE, false, DESC_CUSTOM_SV_FILE);

        SampleVariantsConfig.registerConfig(configBuilder);

        addThreadOptions(configBuilder);

        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        addOutputId(configBuilder);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        addLoggingOptions(configBuilder);
    }
}
