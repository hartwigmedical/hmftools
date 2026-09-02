package com.hartwig.hmftools.geneutils.panelfinder;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.mappability.ProbeQualityProfile.CFG_PROBE_QUALITY_FILE;
import static com.hartwig.hmftools.common.mappability.ProbeQualityProfile.DESC_PROBE_QUALITY_FILE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
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
    public final String GeneIdFile;
    public final RefGenomeVersion RefGenVersion;
    public final String OutputFile;
    public final String OutputBed;

    public final int HighDepthTrimCount;
    public final double MinMappability;
    public final int MinSampleCount;
    public final int GeneUpstreamDistance;
    public final int GeneDownstreamDistance;
    public final boolean RequirePanelGene;

    private static final String HIGH_DEPTH_FILE = "high_depth_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String OUTPUT_BED = "output_bed";

    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String MIN_MAPPABILITY = "min_mappability";
    private static final String MIN_SAMPLE_COUNT = "min_samples";
    private static final String HIGH_DEPTH_TRIM_COUNT = "high_depth_trim_count";
    private static final String GENE_UPSTREAM_DISTANCE = "max_upstream_distance";
    private static final String GENE_DOWNSTREAM_DISTANCE = "max_downstream_distance";
    private static final String REQUIRE_PANEL_GENE = "req_panel_gene";

    protected static final double CHROMOSOME_Y_SAMPLE_FRACTION = 0.4;
    protected static final int DEFAULT_GENE_UPSTREAM_DISTANCE = 10000;
    protected static final int DEFAULT_GENE_DOWNSTREAM_DISTANCE = 0;

    public PanelFinderConfig(final ConfigBuilder configBuilder)
    {
        HighDepthFile = configBuilder.getValue(HIGH_DEPTH_FILE);
        TargetRegionsBed = configBuilder.getValue(TARGET_REGIONS_BED);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        OutputBed = configBuilder.getValue(OUTPUT_BED);
        GeneIdFile = configBuilder.getValue(GENE_ID_FILE);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        EnsemblDataPath = configBuilder.getValue(ENSEMBL_DATA_DIR);
        MappabilityProfileFile = configBuilder.getValue(CFG_PROBE_QUALITY_FILE);

        HighDepthTrimCount = configBuilder.getInteger(HIGH_DEPTH_TRIM_COUNT);
        MinMappability = configBuilder.getDecimal(MIN_MAPPABILITY);
        MinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);
        RequirePanelGene = configBuilder.hasFlag(REQUIRE_PANEL_GENE);

        GeneUpstreamDistance = configBuilder.getInteger(GENE_UPSTREAM_DISTANCE);
        GeneDownstreamDistance = configBuilder.getInteger(GENE_DOWNSTREAM_DISTANCE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(HIGH_DEPTH_FILE, true, "Input regions TSV file - header is optional");
        configBuilder.addPath(TARGET_REGIONS_BED, false, TARGET_REGIONS_BED_DESC);
        configBuilder.addPath(GENE_ID_FILE, false, "Panel gene list");
        EnsemblDataCache.addEnsemblDir(configBuilder, false);
        addRefGenomeVersion(configBuilder);
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output filename");
        configBuilder.addConfigItem(OUTPUT_BED, false, "Output panel definition BED");
        configBuilder.addPath(CFG_PROBE_QUALITY_FILE, false, DESC_PROBE_QUALITY_FILE);

        configBuilder.addDecimal(MIN_MAPPABILITY, "Min mappability to use a high-depth region", 0);
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min high-depth sample count", 0);
        configBuilder.addInteger(HIGH_DEPTH_TRIM_COUNT, "High-depth region trim bases count", 0);
        configBuilder.addInteger(GENE_UPSTREAM_DISTANCE, "Max distance upstream of gene", DEFAULT_GENE_UPSTREAM_DISTANCE);
        configBuilder.addInteger(GENE_DOWNSTREAM_DISTANCE, "Max distance downstream of gene", DEFAULT_GENE_DOWNSTREAM_DISTANCE);
        configBuilder.addFlag(REQUIRE_PANEL_GENE, "Restrict new regions to those overlapping a panel gene");

        addLoggingOptions(configBuilder);
    }

}
