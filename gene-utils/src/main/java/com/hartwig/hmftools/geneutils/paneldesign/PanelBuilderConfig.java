package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.mappability.ProbeQualityProfile.CFG_PROBE_QUALITY_FILE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.Nullable;

public class PanelBuilderConfig
{
    // Reference data files.
    public final String RefGenomeFile;
    public final String EnsemblDir;
    public final String ProbeQualityProfileFile;
    @Nullable
    public final String AmberSitesFile;

    // User input files.
    @Nullable
    public final String TargetGenesFile;
    @Nullable
    public final String CustomRegionsFile;

    // Output config.
    @Nullable
    public final String OutputPrefix;
    public final String OutputDir;
    public final boolean VerboseOutput;

    private static final String CFG_AMBER_SITES_FILE = "amber_sites";
    private static final String DESC_AMBER_SITES_FILE = "Amber het sites file";
    private static final String CFG_TARGET_GENES_FILE = "target_genes";
    private static final String DESC_TARGET_GENES_FILE = "Gene and transcript name file";
    private static final String CFG_CUSTOM_REGIONS_FILE = "custom_regions";
    private static final String DESC_CUSTOM_REGIONS_FILE = "Custom region file";
    private static final String CFG_OUTPUT_PREFIX = "output_prefix";
    private static final String DESC_OUTPUT_PREFIX = "Prefix of output file names";
    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";
    private static final String DESC_VERBOSE_OUTPUT = "Output more information useful for debugging";

    public PanelBuilderConfig(final ConfigBuilder configBuilder)
    {
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        EnsemblDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        ProbeQualityProfileFile = configBuilder.getValue(CFG_PROBE_QUALITY_FILE);
        AmberSitesFile = configBuilder.getValue(CFG_AMBER_SITES_FILE);

        TargetGenesFile = configBuilder.getValue(CFG_TARGET_GENES_FILE);
        CustomRegionsFile = configBuilder.getValue(CFG_CUSTOM_REGIONS_FILE);

        OutputPrefix = configBuilder.hasValue(CFG_OUTPUT_PREFIX) ? configBuilder.getValue(CFG_OUTPUT_PREFIX) : null;
        OutputDir = parseOutputDir(configBuilder);
        VerboseOutput = configBuilder.hasFlag(CFG_VERBOSE_OUTPUT);
    }

    public static void registerConfig(ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        EnsemblDataCache.addEnsemblDir(configBuilder, true);
        ProbeQualityProfile.registerConfig(configBuilder);
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, DESC_AMBER_SITES_FILE);

        configBuilder.addPath(CFG_TARGET_GENES_FILE, false, DESC_TARGET_GENES_FILE);
        configBuilder.addPath(CFG_CUSTOM_REGIONS_FILE, false, DESC_CUSTOM_REGIONS_FILE);

        configBuilder.addConfigItem(CFG_OUTPUT_PREFIX, false, DESC_OUTPUT_PREFIX);
        addOutputDir(configBuilder);
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, DESC_VERBOSE_OUTPUT);

        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
