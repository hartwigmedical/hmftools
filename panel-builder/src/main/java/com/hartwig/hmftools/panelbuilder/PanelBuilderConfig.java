package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.mappability.ProbeQualityProfile.CFG_PROBE_QUALITY_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariantsConfig;

import org.jetbrains.annotations.Nullable;

public record PanelBuilderConfig(
        // Common reference data files.
        String refGenomeFile,
        String ensemblDir,
        String probeQualityProfileFile,
        @Nullable String amberSitesFile,
        // General user input data.
        @Nullable String targetGenesFile,
        @Nullable String customRegionsFile,
        @Nullable SampleVariantsConfig sampleVariants,
        // Output config.
        @Nullable String outputId,
        String outputDir,
        boolean verboseOutput
)
{
    private static final String CFG_AMBER_SITES_FILE = "amber_sites";
    private static final String DESC_AMBER_SITES_FILE = "Amber het sites file";
    private static final String CFG_TARGET_GENES_FILE = "target_genes";
    private static final String DESC_TARGET_GENES_FILE = "Gene and transcript name file";
    private static final String CFG_CUSTOM_REGIONS_FILE = "custom_regions";
    private static final String DESC_CUSTOM_REGIONS_FILE = "Custom region file";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information useful for debugging";

    public static PanelBuilderConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        return new PanelBuilderConfig(
                configBuilder.getValue(REF_GENOME),
                configBuilder.getValue(ENSEMBL_DATA_DIR),
                configBuilder.getValue(CFG_PROBE_QUALITY_FILE),
                configBuilder.getValue(CFG_AMBER_SITES_FILE),
                configBuilder.getValue(CFG_TARGET_GENES_FILE),
                configBuilder.getValue(CFG_CUSTOM_REGIONS_FILE),
                SampleVariantsConfig.fromConfigBuilder(configBuilder),
                configBuilder.getValue(OUTPUT_ID),
                parseOutputDir(configBuilder),
                configBuilder.hasFlag(CFG_VERBOSE_OUTPUT)
        );
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        EnsemblDataCache.addEnsemblDir(configBuilder, true);
        ProbeQualityProfile.registerConfig(configBuilder);
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, DESC_AMBER_SITES_FILE);

        configBuilder.addPath(CFG_TARGET_GENES_FILE, false, DESC_TARGET_GENES_FILE);
        configBuilder.addPath(CFG_CUSTOM_REGIONS_FILE, false, DESC_CUSTOM_REGIONS_FILE);

        addOutputOptions(configBuilder);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        addLoggingOptions(configBuilder);

        SampleVariantsConfig.registerConfig(configBuilder);
    }
}
