package com.hartwig.hmftools.esvee.utils;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_CHR_PARTITION_SIZE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.cli.ParseException;

public class HighDepthConfig
{
    public final String BamFile;
    public final String RefGenome;
    public final String OutputFile;
    public final RefGenomeVersion RefGenVersion;
    public final List<ChrBaseRegion> SpecificRegions;
    public final int PartitionSize;
    public final int HighDepthThreshold;
    public final int Threads;

    private static final String BAM_FILE = "bam_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String HIGH_DEPTH_THRESHOLD = "high_depth_threshold";

    public static final int DEFAULT_HIGH_DEPTH_THRESHOLD = 200;
    public static final int HIGH_DEPTH_REGION_MAX_GAP = 100;

    public HighDepthConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        HighDepthThreshold = configBuilder.getInteger(HIGH_DEPTH_THRESHOLD);

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file to slice for high-depth");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        addRefGenomeConfig(configBuilder, true);

        configBuilder.addInteger(
                HIGH_DEPTH_THRESHOLD, "Level for indicating high-depth", DEFAULT_HIGH_DEPTH_THRESHOLD);

        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
