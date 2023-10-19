package com.hartwig.hmftools.sieve.unmap;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sieve.annotate.AnnotateConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class UnmapperConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(AnnotateConfig.class);

    private static final String BAM_FILE = "bam_file";
    private static final String OUTPUT_BAM_FILE = "output_bam_file";
    private static final String PARTITION_SIZE = "partition_size";

    public final String BamFile;
    public final String OutputBamFile;
    public final int PartitionSize;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    private static final int DEFAULT_PARTITION_SIZE = 100000;

    public UnmapperConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        OutputBamFile = configBuilder.getValue(OUTPUT_BAM_FILE);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file to process");
        configBuilder.addConfigItem(OUTPUT_BAM_FILE, true, "Output BAM file");
        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_PARTITION_SIZE);
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);

        // TODO(m_cooper): Add specific regions config.
        // addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
