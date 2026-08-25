package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.ILLUMINA_DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LOG_READ_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.parseLogReadIds;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SliceParams
{
    public final int TargetDepth;
    public final int MaxPartitionReads;
    public final int MaxRemoteReads;
    public final int ReadLength;
    public final boolean DropDuplicates;
    public final boolean DropConsensus;
    public final List<String> LogReadIds;
    public final boolean CheckLogReads;

    public SliceParams(
            final boolean dropDuplicates, final List<String> logReadIds, final int targetDepth, final int maxPartitionReads,
            final int maxRemoteReads, final int readLength, final boolean dropConsensus)
    {
        TargetDepth = targetDepth;
        MaxPartitionReads = maxPartitionReads;
        MaxRemoteReads = maxRemoteReads;
        ReadLength = readLength;
        DropDuplicates = dropDuplicates;
        DropConsensus = dropConsensus;
        LogReadIds = logReadIds;
        CheckLogReads = !LogReadIds.isEmpty();
    }

    public SliceParams(final ConfigBuilder configBuilder)
    {
        DropDuplicates = configBuilder.hasFlag(DROP_DUPLICATES);
        MaxRemoteReads = configBuilder.getInteger(MAX_REMOTE_READS);
        TargetDepth = configBuilder.getInteger(TARGET_DEPTH);
        MaxPartitionReads = configBuilder.getInteger(MAX_PARTITION_READS);
        ReadLength = configBuilder.getInteger(READ_LENGTH);
        LogReadIds = parseLogReadIds(configBuilder);
        CheckLogReads = !LogReadIds.isEmpty();
        DropConsensus = configBuilder.hasFlag(DROP_CONSENSUS);
    }

    public boolean isDownsampling() { return TargetDepth > 0; }

    private static final String MAX_PARTITION_READS = "max_partition_reads";
    private static final String MAX_REMOTE_READS = "max_remote_reads";
    private static final String TARGET_DEPTH = "target_depth";
    private static final String DROP_DUPLICATES = "drop_duplicates";
    private static final String READ_LENGTH = "read_length";
    private static final String DROP_CONSENSUS = "drop_consensus";

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(MAX_PARTITION_READS, "Max partition reads (perf-only)", 0);
        configBuilder.addInteger(MAX_REMOTE_READS, "Max remote reads (perf-only)", 0);
        configBuilder.addInteger(TARGET_DEPTH, "Target depth impacting down-sampling (zero means not applied)", 0);
        configBuilder.addFlag(DROP_DUPLICATES, "Ignore duplicate reads");
        configBuilder.addInteger(READ_LENGTH,  "Read length for use in down-sampled target depth", ILLUMINA_DEFAULT_READ_LENGTH);
        configBuilder.addConfigItem(LOG_READ_IDS, LOG_READ_IDS_DESC);
        configBuilder.addFlag(DROP_CONSENSUS, "Drop consensus reads from Redux input BAM");
    }
}
