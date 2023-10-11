package com.hartwig.hmftools.sieve.bedstats;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BedStatsConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(BedStatsConfig.class);

    private static final String BAM_FILE = "bam_file";
    private static final String HIGH_DEPTH_FILE = "high_depth_file";
    private static final String OUTPUT_FILE = "output_file";

    public final String BamFile;
    public final String HighDepthFile;
    public final String OutputFile;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    public BedStatsConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        HighDepthFile = configBuilder.getValue(HIGH_DEPTH_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file");
        configBuilder.addPath(HIGH_DEPTH_FILE, true, "high depth file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
    }
}
