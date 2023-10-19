package com.hartwig.hmftools.sieve.unmap;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

// TODO(m_cooper): Unit tests.
// TODO(m_cooper): Write readme.
public class Unmapper
{
    private final UnmapperConfig mConfig;

    public Unmapper(final ConfigBuilder configBuilder)
    {
        mConfig = new UnmapperConfig(configBuilder);
    }

    // TODO(m_cooper): Use ConcurrentHashMap?
    public static synchronized void registerReadForDropping(final Set<PrimaryReadInfo> readsToDrop, final PrimaryReadInfo read)
    {
        readsToDrop.add(read);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        UnmapperConfig.addConfig(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        // TODO(m_cooper): Fix null pointer exception here.
        //        logVersion();

        Unmapper mapDropper = new Unmapper(configBuilder);
        mapDropper.run();
    }

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("mark-dups.version");
        MD_LOGGER.info("MarkDups version: {}", version.version());
    }

    public void run()
    {
        if(mConfig.BamFile == null)
        {
            MD_LOGGER.error("no BAM file specified");
            System.exit(1);
        }

        if(mConfig.OutputBamFile == null)
        {
            MD_LOGGER.error("no output BAM file specified");
            System.exit(1);
        }

        // TODO(m_cooper): Use ConcurrentLinkedQueue?
        final List<ChrBaseRegion> allPartitions = Lists.newArrayList();
        for(final HumanChromosome chromosome : HumanChromosome.values())
        {
            final String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
            allPartitions.addAll(partitionChromosome(chromosomeStr, mConfig.RefGenVersion, mConfig.PartitionSize));
        }

        final ArrayBlockingQueue<ChrBaseRegion> jobs = new ArrayBlockingQueue<>(allPartitions.size(), true, allPartitions);
        final Set<PrimaryReadInfo> readsToDrop = Sets.newHashSet();
        final List<UnmapCollector> collectorTasks = Lists.newArrayList();
        for(int i = 0; i < Math.max(1, mConfig.Threads); ++i)
        {
            collectorTasks.add(new UnmapCollector(mConfig, jobs, readsToDrop));
        }

        final List<Callable> callableList = collectorTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        MD_LOGGER.info("unmapper complete");
    }
}
