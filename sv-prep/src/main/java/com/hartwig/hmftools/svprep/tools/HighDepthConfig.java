package com.hartwig.hmftools.svprep.tools;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_CHR_PARTITION_SIZE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
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
    private static final String THREADS = "threads";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String HIGH_DEPTH_THRESHOLD = "high_depth_threshold";

    private static final int DEFAULT_HIGH_DEPTH_THRESHOLD = 200;

    public HighDepthConfig(final CommandLine cmd)
    {
        BamFile = cmd.getOptionValue(BAM_FILE);
        OutputFile = cmd.getOptionValue(OUTPUT_FILE);
        RefGenome = cmd.getOptionValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));
        HighDepthThreshold = Integer.parseInt(cmd.getOptionValue(HIGH_DEPTH_THRESHOLD, String.valueOf(DEFAULT_HIGH_DEPTH_THRESHOLD)));

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }
    }

    public static void addOptions(final Options options)
    {
        options.addOption(BAM_FILE, true, "BAM file to slice for high-depth");
        options.addOption(OUTPUT_FILE, true, "Output file");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(HIGH_DEPTH_THRESHOLD, true, "Level for indicating high-depth");
        options.addOption(PARTITION_SIZE, true, "Partition size, default = " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(THREADS, true, "Multi-thread count");
        addSpecificChromosomesRegionsConfig(options);
    }
}
