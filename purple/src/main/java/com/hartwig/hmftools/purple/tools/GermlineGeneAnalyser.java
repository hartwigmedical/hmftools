package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.perf.TaskExecutor;

import org.jetbrains.annotations.NotNull;

public class GermlineGeneAnalyser
{
    static final String PURPLE_GERMLINE_GENE_DATA_CSV = "purple_germline_gene_data.csv";
    static final String MIN_FREQUENCY_CFG = "min_frequency";
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final int mThreads;
    private final int mMinFrequency;
    private final String mOutputDir;
    private final RefGenomeVersion mRefGenomeVersion;

    public GermlineGeneAnalyser(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDataDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);
        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);

        mOutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        mMinFrequency = configBuilder.getInteger(MIN_FREQUENCY_CFG);
    }

    public void runAnalysis()
    {
        removeInvalidSampleIds();
        PPL_LOGGER.info("running Purple germline gene analysis for {} samples", mSampleIds.size());

        List<SampleGermlineGeneTask> sampleTasks = createAnalysisTasks();
        PPL_LOGGER.info("tasks created: {}", sampleTasks.size());

        final List<Callable<Void>> callableList = new ArrayList<>(sampleTasks);
        TaskExecutor.executeTasks(callableList, mThreads);
        PPL_LOGGER.info("files read");
        EventCounts results = new EventCounts(new TreeMap<>());
        sampleTasks.forEach(sampleTask -> results.mergeAdd(sampleTask.getEventCounts()));

        String fileName = mOutputDir + "/" + PURPLE_GERMLINE_GENE_DATA_CSV;
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);
            results.writeTo(writer, mMinFrequency, mRefGenomeVersion);
            writer.close();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("could not write output to " + fileName, e);
        }
        PPL_LOGGER.info("Purple germline event frequency analysis complete");
    }

    @NotNull
    private List<SampleGermlineGeneTask> createAnalysisTasks()
    {
        List<SampleGermlineGeneTask> sampleTasks = Lists.newArrayList();
        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleGermlineGeneTask(i, mPurpleDataDir));
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            if(taskIndex >= sampleTasks.size())
            {
                taskIndex = 0;
            }
            sampleTasks.get(taskIndex).getSampleIds().add(sampleId);
            ++taskIndex;
        }
        return sampleTasks;
    }

    private void removeInvalidSampleIds()
    {
        if(mSampleIds.isEmpty())
        {
            PPL_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        List<String> invalidSampleIds = findInvalidSampleIds();
        if(!invalidSampleIds.isEmpty())
        {
            PPL_LOGGER.warn("invalid sampleIds: {}", invalidSampleIds);
        }
        mSampleIds.removeAll(invalidSampleIds);
    }

    private List<String> findInvalidSampleIds()
    {
        List<String> result = new ArrayList<>();
        for(String sampleId : mSampleIds)
        {
            File segmentFile = new File(GermlineAmpDel.generateFilename(mPurpleDataDir, sampleId));
            if(!segmentFile.exists())
            {
                result.add(sampleId);
            }
        }
        return result;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, true);
        configBuilder.addPath(PURPLE_DIR_CFG, true, "Directory pattern for sample purple directory");
        addLoggingOptions(configBuilder);
        addRefGenomeVersion(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
        configBuilder.addInteger(MIN_FREQUENCY_CFG, "lowest frequency of germline events to include in output", 4);
        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);
        GermlineGeneAnalyser germlineGeneAnalyser = new GermlineGeneAnalyser(configBuilder);
        germlineGeneAnalyser.runAnalysis();
    }
}