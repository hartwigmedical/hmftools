package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.purple.HrdData;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class HrdDetectionAnalyser
{
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final String mChordDir;
    private final int mThreads;

    private final BufferedWriter mWriter;

    public HrdDetectionAnalyser(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
            mSampleIds = Lists.newArrayList(configBuilder.getValue(SAMPLE));
        else
            mSampleIds = loadSampleIdsFile(configBuilder);

        mPurpleDataDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mChordDir = configBuilder.getValue(CHORD_DIR_CFG);
        mThreads = configBuilder.getInteger(THREADS);

        String outputId = configBuilder.getValue(OUTPUT_ID);
        String outputDir = parseOutputDir(configBuilder);
        mWriter = initialiseWriter(outputDir, outputId);
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            PPL_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        if(mPurpleDataDir == null)
        {
            PPL_LOGGER.error("missing purple & chord data paths");
            System.exit(1);
        }

        PPL_LOGGER.info("running Purple HRD analysis for {} samples", mSampleIds.size());

        List<SampleTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleTask(i));
        }

        // allocate samples to threads
        int taskIndex = 0;
        for(String sample : mSampleIds)
        {
            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;

            sampleTasks.get(taskIndex).SampleIds.add(sample);

            ++taskIndex;
        }

        final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        closeBufferedWriter(mWriter);

        PPL_LOGGER.info("Purple HRD analysis complete");
    }

    private class SampleTask implements Callable<Void>
    {
        private final int mTaskId;
        public final List<String> SampleIds;

        public SampleTask(int taskId)
        {
            mTaskId = taskId;
            SampleIds = Lists.newArrayList();
        }

        @Override
        public Void call()
        {
            for(int i = 0; i < SampleIds.size(); ++i)
            {
                processSample(SampleIds.get(i));

                if((i % 100) == 0)
                    PPL_LOGGER.info("{}: processed {} samples", mTaskId, i);

            }

            return null;
        }

        private void processSample(final String sampleId)
        {
            List<PurpleCopyNumber> copyNumbers = null;
            PurityContext purityContext = null;
            ChordData chordData = null;

            try
            {
                if(mChordDir != null)
                {
                    String sampleChordDirectory = convertWildcardSamplePath(mChordDir, sampleId);
                    chordData = ChordDataFile.read(ChordDataFile.generateFilename(sampleChordDirectory, sampleId));
                }
                else
                {
                    chordData = ImmutableChordData.builder()
                            .BRCA1Value(0)
                            .BRCA2Value(0)
                            .hrdType("N/A")
                            .hrdValue(0)
                            .hrStatus(ChordStatus.UNKNOWN)
                            .remarksHrdType("")
                            .remarksHrStatus("").build();
                }

                String samplePurpleDirectory = convertWildcardSamplePath(mPurpleDataDir, sampleId);
                copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(samplePurpleDirectory, sampleId));
                purityContext = PurityContextFile.read(samplePurpleDirectory, sampleId);
            }
            catch(IOException e)
            {
                PPL_LOGGER.error("failed to load file data {}", e.toString());
                return;
            }

            if(chordData == null || purityContext == null)
            {
                PPL_LOGGER.info("sample({}) invalid purity({}) or chord({}) data",
                        sampleId, purityContext != null ? "valid" : "missing", chordData != null ? "valid" : "missing");
                return;
            }

            PPL_LOGGER.debug(format("sample(%s) cnRecords(%d) chord(%s %.3f)",
                    sampleId, copyNumbers.size(), chordData.hrStatus(), chordData.hrdValue()));

            HrdDetection hrdDetection = new HrdDetection();
            final HrdData hrdData = hrdDetection.calculateHrdData(copyNumbers);

            writeSampleData(sampleId, hrdData, purityContext, chordData);
        }
    }

    private synchronized BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String fileName = outputDir;

            if(mSampleIds.size() == 1)
            {
                fileName += mSampleIds.get(0) + ".purple.hrd.tsv";
            }
            else
            {
                fileName += "cohort_purple_hrd_analysis";

                if(outputId != null)
                    fileName += "." + outputId;

                fileName += TSV_EXTENSION;
            }

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mSampleIds.size() > 1)
                sj.add("SampleId");

            sj.add("Purity").add("Ploidy");
            sj.add("LohSegments").add("SegmentBreaks").add("SegmentImbalances").add("HrdStatus");

            if(mChordDir != null)
            {
                sj.add("ChordHrdStatus").add("ChordHrd").add("ChordBRCA1").add("ChordBRCA2");
            }

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeSampleData(
            final String sampleId, final HrdData hrdData, final PurityContext purityContext, final ChordData chordData)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mSampleIds.size() > 1)
            {
                sj.add(sampleId);
            }

            sj.add(String.valueOf(purityContext.bestFit().purity()));
            sj.add(String.valueOf(purityContext.bestFit().ploidy()));

            sj.add(String.valueOf(hrdData.LohSegments));
            sj.add(String.valueOf(hrdData.SegmentBreaks));
            sj.add(String.valueOf(hrdData.SegmentImbalances));
            sj.add(hrdData.Status.toString());

            if(mChordDir != null)
            {
                sj.add(chordData.hrStatus().toString());
                sj.add(format("%.3f", chordData.hrdValue()));
                sj.add(format("%.3f", chordData.BRCA1Value()));
                sj.add(format("%.3f", chordData.BRCA2Value()));
            }

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene overlap file output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(CHORD_DIR_CFG, false, CHORD_DIR_DESC);

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HrdDetectionAnalyser hrdDetectionAnalyser = new HrdDetectionAnalyser(configBuilder);
        hrdDetectionAnalyser.run();
    }
}
