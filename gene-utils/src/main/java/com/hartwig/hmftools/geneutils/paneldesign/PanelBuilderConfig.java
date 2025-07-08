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
    public final String GeneTranscriptFile;
    public final String CustomRegionFile;

    public final String RefGenomeFile;
    public final String EnsemblDir;
    public final String ProbeQualityProfileFile;

    public final String OutputPrefix;
    public final String OutputDir;

    private static final String CFG_AMBER_SITES_FILE = "amber_sites_file";
    private static final String CFG_GENE_TRANSCRIPT_FILE = "gene_transcript_file";
    private static final String CFG_CUSTOM_REGION_FILE = "custom_region_file";
    private static final String CFG_OUTPUT_PREFIX = "output_prefix";

    public PanelBuilderConfig(final ConfigBuilder configBuilder)
    {
        AmberSitesFile = configBuilder.getValue(CFG_AMBER_SITES_FILE);
        GeneTranscriptFile = configBuilder.getValue(CFG_GENE_TRANSCRIPT_FILE);
        CustomRegionFile = configBuilder.getValue(CFG_CUSTOM_REGION_FILE);

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
        configBuilder.addPath(CFG_AMBER_SITES_FILE, false, "Amber het sites file");
        configBuilder.addPath(CFG_GENE_TRANSCRIPT_FILE, false, "Gene and transcript name file");
        configBuilder.addPath(CFG_CUSTOM_REGION_FILE, false, "Custom region file");

        EnsemblDataCache.addEnsemblDir(configBuilder, true);

        addRefGenomeFile(configBuilder, true);

        ProbeQualityProfile.registerConfig(configBuilder);

        configBuilder.addConfigItem(CFG_OUTPUT_PREFIX, true, "prefix of output BED filename");
        addOutputDir(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
