package com.hartwig.hmftools.sieve;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sieve.MapDropperConfig.MD_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

// TODO(m_cooper): Unit tests.
// TODO(m_cooper): Write readme.
public class MapDropper
{
    private final MapDropperConfig mConfig;

    public MapDropper(@NotNull final ConfigBuilder configBuilder)
    {
        mConfig = new MapDropperConfig(configBuilder);
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

        final RecordWriter writer = new RecordWriter(mConfig);

        // TODO(m_cooper): Multi-thread.
        final List<MapDropperTask> dropTasks = Lists.newArrayList();
        dropTasks.add(new MapDropperTask(mConfig, writer));
        final List<Callable> callableList = dropTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        writer.close();

        MD_LOGGER.info("map dropper complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        MapDropperConfig.addConfig(configBuilder);
        addLoggingOptions(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        // TODO(m_cooper): Fix null pointer exception here.
//        logVersion();

        MapDropper mapDropper = new MapDropper(configBuilder);
        mapDropper.run();
    }

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("mark-dups.version");
        MD_LOGGER.info("MarkDups version: {}", version.version());
    }
}
