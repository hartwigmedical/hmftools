package com.hartwig.hmftools.bamtools.bamtofastq;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// TODO NEXT: TEST
public class BamToFastqConfig
{
    public static final Logger BFQ_LOGGER = LogManager.getLogger(BamToFastqConfig.class);

    private static final String BAM_FILE = "bam_file";
    public static final String FASTQ1_OUTPUT_FILE = "fastq1";
    public static final String FASTQ2_OUTPUT_FILE = "fastq2";
    public static final String NO_WRITE = "no_write";
    private static final String KEEP_CONSENSUS_READS = "keep_consensus_reads";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String SILENT_VALIDATION = "silent_validation";
    private static final String PERF_DEBUG = "perf_debug";

    private static final int DEFAULT_PARTITION_SIZE = 1_000_000;

    public final String RefGenomeFile;
    public final RefGenomeSource RefGenome;
    public final String BamFile;
    public final String Fastq1OutputFile;
    public final String Fastq2OutputFile;
    public final int Threads;
    public final boolean NoWrite;
    public final boolean KeepConsensusReads;
    public final int PartitionSize;
    public final boolean SilentValidation;
    public final boolean PerfDebug;

    public BamToFastqConfig(final ConfigBuilder configBuilder)
    {
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);
        BamFile = configBuilder.getValue(BAM_FILE);
        Fastq1OutputFile = configBuilder.getValue(FASTQ1_OUTPUT_FILE);
        Fastq2OutputFile = configBuilder.getValue(FASTQ2_OUTPUT_FILE);
        Threads = parseThreads(configBuilder);
        NoWrite = configBuilder.hasFlag(NO_WRITE);
        KeepConsensusReads = configBuilder.hasFlag(KEEP_CONSENSUS_READS);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        SilentValidation = configBuilder.hasFlag(SILENT_VALIDATION);
        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);

        validate();
    }

    public boolean hasFastqOutputFiles()
    {
        return Fastq1OutputFile != null && Fastq2OutputFile != null;
    }

    private void validate()
    {
        if(!(NoWrite || hasFastqOutputFiles()))
        {
            BFQ_LOGGER.error("Must specify at either -{} or -{} and -{}", NO_WRITE, FASTQ1_OUTPUT_FILE, FASTQ2_OUTPUT_FILE);
            System.exit(1);
        }

        if(NoWrite && hasFastqOutputFiles())
        {
            BFQ_LOGGER.error("If specifying -{} then cannot specify -{} -{}", NO_WRITE, FASTQ1_OUTPUT_FILE, FASTQ2_OUTPUT_FILE);
            System.exit(1);
        }
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC);
        configBuilder.addPath(BAM_FILE, true, "Input BAM");
        configBuilder.addConfigItem(FASTQ1_OUTPUT_FILE, false, "Output path for first in pair fastq entries");
        configBuilder.addConfigItem(FASTQ2_OUTPUT_FILE, false, "Output path for second in pair fastq entries");
        configBuilder.addFlag(NO_WRITE, "Do not write output");
        configBuilder.addFlag(KEEP_CONSENSUS_READS, "Keep consensus reads");
        configBuilder.addInteger(PARTITION_SIZE, "Size of chromosome partitions that are consumed by threaded workers", DEFAULT_PARTITION_SIZE);
        configBuilder.addFlag(SILENT_VALIDATION, "Do not validate SAM records upon reading");
        configBuilder.addFlag(PERF_DEBUG, "Logs performance statistics");

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    @VisibleForTesting
    public BamToFastqConfig(boolean keepConsensusReads)
    {
        RefGenomeFile = null;
        RefGenome = null;
        BamFile = null;
        Fastq1OutputFile = null;
        Fastq2OutputFile = null;
        Threads = 1;
        NoWrite = false;
        KeepConsensusReads = keepConsensusReads;
        PartitionSize = DEFAULT_PARTITION_SIZE;
        SilentValidation = false;
        PerfDebug = false;
    }
}
