package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.APP_NAME;
import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.common.MutContextCount;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class ChordDataPrep
{
    public final ChordConfig mConfig;

    public ChordDataPrep(final ConfigBuilder configBuilder) { mConfig = new ChordConfig(configBuilder); }

    public ChordDataPrep(final ChordConfig config) { mConfig = config; }

    public void prepSingleSample() throws IOException
    {
        String sampleId = mConfig.SampleIds.get(0);
        CHORD_LOGGER.info("Extracting CHORD features for sample: {}", sampleId);

        SamplePrepTask sampleTask = new SamplePrepTask(mConfig, 0, null);
        sampleTask.processSample();

        List<MutContextCount> contextCounts = sampleTask.mContextCounts;

        ChordDataWriter writer = new ChordDataWriter(mConfig);
        writer.writeHeader(contextCounts);
        writer.writeValues(sampleId, contextCounts);
        writer.close();
    }

    public void prepMultiSample() throws IOException
    {
        CHORD_LOGGER.info("Extracting CHORD features in multi sample mode: {} samples, {} threads", mConfig.SampleIds.size(), mConfig.Threads);

        ConcurrentHashMap<String, List<MutContextCount>> contextCountsMatrix = new ConcurrentHashMap<>();

        List<SamplePrepTask> sampleTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); sampleIndex++)
        {
            SamplePrepTask sampleTask = new SamplePrepTask(mConfig, sampleIndex, contextCountsMatrix);
            sampleTasks.add(sampleTask);
        }

        List<Callable> callableTasks = new ArrayList<>(sampleTasks);
        TaskExecutor.executeTasks(callableTasks, mConfig.Threads);

        ChordDataWriter writer = new ChordDataWriter(mConfig);
        writer.writeHeader(contextCountsMatrix.get(mConfig.SampleIds.get(0)));

        for(String sampleId : mConfig.SampleIds)
        {
            List<MutContextCount> sampleContextCounts = contextCountsMatrix.get(sampleId);
            writer.writeValues(sampleId, sampleContextCounts);
        }

        writer.close();
    }

    public void run()
    {
        if(mConfig.SampleIds.isEmpty())
        {
            CHORD_LOGGER.error("No sample ID(s) loaded");
            System.exit(1);
        }

        try
        {
            long startTimeMs = System.currentTimeMillis();

            if(mConfig.isSingleSample())
                prepSingleSample();
            else
                prepMultiSample();

            CHORD_LOGGER.info("CHORD feature extraction complete, mins({})", runTimeMinsStr(startTimeMs));
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("CHORD feature extraction failed: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ChordConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ChordDataPrep chordDataPrep = new ChordDataPrep(configBuilder);
        chordDataPrep.run();
    }
}
