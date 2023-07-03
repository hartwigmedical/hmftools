package com.hartwig.hmftools.purple.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HrdDetectionAnalyser
{
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final String mChordDir;
    private final int mThreads;

    private final HrdDetection mHrdDetection;

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

        mHrdDetection = new HrdDetection();

        mWriter = initialiseWriter(parseOutputDir(configBuilder));
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

        int processed = 0;

        for(String sampleId : mSampleIds)
        {
            processSample(sampleId);

            ++processed;

            if((processed % 100) == 0)
                PPL_LOGGER.info("processed {} samples", processed);
        }

        closeBufferedWriter(mWriter);

        PPL_LOGGER.info("Purple HRD analysis complete");
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
                String sampleChordDirectory =  convertWildcardSamplePath(mChordDir, sampleId);
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

        PPL_LOGGER.debug(format("sample(%s) ploidy(%.1f) cnRecords(%d) chord(%s %.3f)",
                sampleId, purityContext.bestFit().ploidy(), copyNumbers.size(), chordData.hrStatus(), chordData.hrdValue()));

        final HrdData hrdData = mHrdDetection.calculateHrdData(copyNumbers, purityContext.bestFit().ploidy());

        writeSampleData(sampleId, chordData, hrdData);
    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String fileName = outputDir;

            if(mSampleIds.size() == 1)
                fileName += mSampleIds.get(0) + ".purple.hrd.tsv";
            else
                fileName += "cohort_purple_hrd_analysis.tsv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mSampleIds.size() == 1)
                writer.write("SampleId\t");

            writer.write("HRDStatus\tHRD\tBRCA1\tBRCA2");
            writer.write("\tLohSegments\tSegmentBreaks\tSegmentImbalances\tScore");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeSampleData(final String sampleId, final ChordData chordData, final HrdData hrdData)
    {
        try
        {
            if(mSampleIds.size() == 1)
            {
                mWriter.write(format("%s\t", sampleId));
            }

            mWriter.write(format("%s\t%.3f\t%.3f\t%.3f",
                    chordData.hrStatus(), chordData.hrdValue(), chordData.BRCA1Value(), chordData.BRCA2Value()));

            mWriter.write(format("\t%d\t%.1f\t%d\t%.1f",
                    hrdData.LohSegments, hrdData.SegmentBreaks, hrdData.SegmentImbalances, hrdData.score()));

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

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        HrdDetectionAnalyser hrdDetectionAnalyser = new HrdDetectionAnalyser(configBuilder);
        hrdDetectionAnalyser.run();
    }
}
