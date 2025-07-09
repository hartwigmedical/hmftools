package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.nio.file.Paths;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class PanelBuilderConfig
{
    public final String AmberSitesFile;
    public final String TargetGenesFile;
    public final String CustomRegionsFile;

    public final String RefGenomeFile;
    public final String EnsemblDir;
    public final String ProbeQualityProfileFile;

    public final String OutputPrefix;
    public final String OutputDir;

    private static final String CFG_AMBER_SITES_FILE = "amber_sites_file";
    private static final String DESC_AMBER_SITES_FILE = "Amber het sites file";
    private static final String CFG_TARGET_GENES_FILE = "gene_transcript_file";
    private static final String DESC_TARGET_GENES_FILE = "Gene and transcript name file";
    private static final String CFG_CUSTOM_REGIONS_FILE = "custom_region_file";
    private static final String DESC_CUSTOM_REGIONS_FILE = "Custom region file";
    private static final String CFG_OUTPUT_PREFIX = "output_prefix";
    private static final String DESC_OUTPUT_PREFIX = "Prefix of output file names";

    public PanelBuilderConfig(final ConfigBuilder configBuilder)
    {
        AmberSitesFile = configBuilder.getValue(CFG_AMBER_SITES_FILE);
        TargetGenesFile = configBuilder.getValue(CFG_TARGET_GENES_FILE);
        CustomRegionsFile = configBuilder.getValue(CFG_CUSTOM_REGIONS_FILE);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        EnsemblDir = configBuilder.getValue(EnsemblDataCache.ENSEMBL_DATA_DIR);

        ProbeQualityProfileFile = configBuilder.getValue(ProbeQualityProfile.PROBE_QUALITY_FILE_CONFIG);

        OutputPrefix = configBuilder.getValue(CFG_OUTPUT_PREFIX);
        OutputDir = parseOutputDir(configBuilder);
    }

    public String outputFilePath(final String fileName)
    {
        return Paths.get(OutputDir, format("%s.%s", OutputPrefix, fileName)).toString();
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, DESC_AMBER_SITES_FILE);
        configBuilder.addPath(CFG_TARGET_GENES_FILE, false, DESC_TARGET_GENES_FILE);
        configBuilder.addPath(CFG_CUSTOM_REGIONS_FILE, false, DESC_CUSTOM_REGIONS_FILE);

        EnsemblDataCache.addEnsemblDir(configBuilder, true);

        addRefGenomeFile(configBuilder, true);

        ProbeQualityProfile.registerConfig(configBuilder);

        configBuilder.addConfigItem(CFG_OUTPUT_PREFIX, true, DESC_OUTPUT_PREFIX);
        addOutputDir(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
