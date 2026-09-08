package com.hartwig.hmftools.bamtools.depth;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.GenicRegions.FIXED_GENE_REGIONS;
import static com.hartwig.hmftools.bamtools.depth.GenicRegions.REMOVE_GENE_OVERLAPS;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.commons.cli.ParseException;

public class CombinerConfig
{
    public final List<String> SampleIds;
    public final String HighDepthFiles;

    public final int Threads;
    public final int MinSampleCount;
    public final int MinRegionSize;
    public final int ReadLength;
    public final RefGenomeVersion RefGenVersion;

    public final List<ChrBaseRegion> SpecificRegions;
    public final String OutputFile;
    public final boolean WriteWithLabel;

    // config
    private static final String HIGH_DEPTH_FILES = "high_depth_files";
    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";
    private static final String MIN_REGION_SIZE = "min_region_size";
    private static final String WRITE_LABEL = "write_label";
    private static final String READ_LENGTH = "read_length";

    // constants
    private static final int DEFAULT_MIN_SAMPLE_COUNT = 4;

    public CombinerConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = loadSampleIdsFile(configBuilder);
        HighDepthFiles = configBuilder.getValue(HIGH_DEPTH_FILES);

        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        MinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);
        MinRegionSize = configBuilder.getInteger(MIN_REGION_SIZE);
        ReadLength = configBuilder.getInteger(READ_LENGTH);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        WriteWithLabel = configBuilder.hasFlag(WRITE_LABEL);

        Threads = parseThreads(configBuilder);

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            BT_LOGGER.error("failed to load specific regions");
            System.exit(1);
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(HIGH_DEPTH_FILES, true, "High depth sample file(s), use '*' in for sampleId");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        configBuilder.addPath(FIXED_GENE_REGIONS, false, "Reference blacklist file to include");
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min sample count to produce region", DEFAULT_MIN_SAMPLE_COUNT);
        configBuilder.addInteger(MIN_REGION_SIZE, "Min final region width", 0);
        configBuilder.addInteger(READ_LENGTH, "BAM read length", 0);
        configBuilder.addFlag(REMOVE_GENE_OVERLAPS, "Remove high depth regions that overlap driver or fusion genes");
        configBuilder.addFlag(WRITE_LABEL, "Write depth info as 'Label' column for compatibility with panel definition");
        configBuilder.addConfigItem(REF_GENOME_VERSION, REF_GENOME_VERSION_CFG_DESC);

        addGenePanelOption(configBuilder, false);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
