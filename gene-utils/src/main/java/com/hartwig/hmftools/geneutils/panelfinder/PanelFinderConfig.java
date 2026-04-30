package com.hartwig.hmftools.geneutils.panelfinder;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PanelFinderConfig
{
    public final String HighDepthFile;
    public final String TargetRegionsBed;
    public final String EnsemblDataPath;
    public final String MappabilityProfileFile;
    public final RefGenomeVersion RefGenVersion;
    public final String OutputFile;
    public final String OutputBed;

    public final int HighDepthTrimCount;
    public final double MinMappability;
    public final int MinSampleCount;

    private static final String HIGH_DEPTH_FILE = "high_depth_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String OUTPUT_BED = "output_bed";

    private static final String MIN_MAPPABILITY = "min_mappability";
    private static final String MIN_SAMPLE_COUNT = "min_samples";
    private static final String HIGH_DEPTH_TRIM_COUNT = "high_depth_trim_count";

    protected static final double CHROMOSOME_Y_SAMPLE_FRACTION = 0.4;

    public PanelFinderConfig(final ConfigBuilder configBuilder)
    {
        HighDepthFile = configBuilder.getValue(HIGH_DEPTH_FILE);
        TargetRegionsBed = configBuilder.getValue(TARGET_REGIONS_BED);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        OutputBed = configBuilder.getValue(OUTPUT_BED);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        EnsemblDataPath = configBuilder.getValue(ENSEMBL_DATA_DIR);
        MappabilityProfileFile = configBuilder.getValue(ProbeQualityProfile.CFG_PROBE_QUALITY_FILE);

        HighDepthTrimCount = configBuilder.getInteger(HIGH_DEPTH_TRIM_COUNT);
        MinMappability = configBuilder.getDecimal(MIN_MAPPABILITY);
        MinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(HIGH_DEPTH_FILE, true, "Input regions TSV file - header is optional");
        configBuilder.addPath(TARGET_REGIONS_BED, false, TARGET_REGIONS_BED_DESC);
        EnsemblDataCache.addEnsemblDir(configBuilder, false);
        addRefGenomeVersion(configBuilder);
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output filename");
        configBuilder.addConfigItem(OUTPUT_BED, false, "Output panel definition BED");
        ProbeQualityProfile.registerConfig(configBuilder);

        configBuilder.addDecimal(MIN_MAPPABILITY, "Min mappability to use a high-depth region", 0);
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min high-depth sample count", 0);
        configBuilder.addInteger(HIGH_DEPTH_TRIM_COUNT, "High-depth region trim bases count", 0);

        addLoggingOptions(configBuilder);
    }

}
