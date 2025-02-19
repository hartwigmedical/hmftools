package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class JitterAnalyserConfig
{
    public final String SampleId;
    public final RefGenomeVersion RefGenVersion;
    public final String RefGenomeFile;
    public final SequencingType Sequencing;

    public final String RefGenomeMsiFile;
    public final String OutputDir;

    public final int MinMappingQuality;
    public final int MaxSitesPerType;
    public final double MaxSingleSiteAltContribution;

    public final boolean WritePlots;
    public final boolean WriteSiteFile;

    public static final String JITTER_WRITE_MSI_PLOTS = "write_msi_plots";
    public static final String JITTER_WRITE_MSI_PLOTS_DESC = "Generate MSI jitter plots";

    public static final String JITTER_MSI_SITES_FILE = "ref_genome_msi_file";
    public static final String JITTER_MSI_SITES_FILE_DESC = "Path to ref genome MSI sites file";

    public static final String JITTER_MAX_SITES_PER_TYPE = "msi_max_sites_per_type";
    public static final String JITTER_MAX_SITES_PER_TYPE_DESC = "Max number of sites per microsatellite unit / length type";

    private static final String MIN_MAP_QUALITY = "msi_min_map_quality";
    private static final String MAX_SINGLE_SITE_ALT_CONTRIBUTION = "max_site_alt_contribution";

    public static final String JITTER_WRITE_SITE_FILE = "write_msi_site_file";
    public static final String JITTER_WRITE_SITE_FILE_DESC = "Write MSI site file (useful for debugging)";

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_NUM_SITES_PER_TYPE = 5_000;

    private JitterAnalyserConfig(final String sampleId, final String refGenomeFile, final RefGenomeVersion refGenVersion,
            final SequencingType sequencing, final String outputDir, final ConfigBuilder configBuilder)
    {
        SampleId = sampleId;
        RefGenomeFile = refGenomeFile;
        RefGenVersion = refGenVersion;
        Sequencing = sequencing;
        OutputDir = outputDir;

        RefGenomeMsiFile = configBuilder.getValue(JITTER_MSI_SITES_FILE);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(JITTER_MAX_SITES_PER_TYPE);
        MaxSingleSiteAltContribution = configBuilder.getDecimal(MAX_SINGLE_SITE_ALT_CONTRIBUTION);
        WritePlots = configBuilder.hasFlag(JITTER_WRITE_MSI_PLOTS);
        WriteSiteFile = configBuilder.hasFlag(JITTER_WRITE_SITE_FILE);
    }

    @Nullable
    public static JitterAnalyserConfig create(final String sampleId, final String refGenomeFile, final RefGenomeVersion refGenVersion,
            final SequencingType sequencing, final String outputDir, final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(JITTER_MSI_SITES_FILE))
            return null;

        return new JitterAnalyserConfig(sampleId, refGenomeFile, refGenVersion, sequencing, outputDir, configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(JITTER_MSI_SITES_FILE, false, JITTER_MSI_SITES_FILE_DESC);
        configBuilder.addInteger(MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);
        configBuilder.addInteger(JITTER_MAX_SITES_PER_TYPE, JITTER_MAX_SITES_PER_TYPE_DESC, DEFAULT_NUM_SITES_PER_TYPE);
        configBuilder.addDecimal(MAX_SINGLE_SITE_ALT_CONTRIBUTION, "Max percentage a single alt site can contribute", DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION);
        configBuilder.addFlag(JITTER_WRITE_MSI_PLOTS, JITTER_WRITE_MSI_PLOTS_DESC);
        configBuilder.addFlag(JITTER_WRITE_SITE_FILE, JITTER_WRITE_SITE_FILE_DESC);
    }
}
