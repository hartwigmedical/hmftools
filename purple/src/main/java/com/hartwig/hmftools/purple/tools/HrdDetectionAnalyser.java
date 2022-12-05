package com.hartwig.hmftools.purple.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;

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

    private static final String PURPLE_DATA_DIR = "purple_dir";
    private static final String CHORD_DIR = "chord_dir";

    public HrdDetectionAnalyser(final CommandLine cmd)
    {
        mSampleIds = loadSampleIdsFile(cmd);
        mPurpleDataDir = cmd.getOptionValue(PURPLE_DATA_DIR);
        mChordDir = cmd.getOptionValue(CHORD_DIR);
        mThreads = parseThreads(cmd);

        mHrdDetection = new HrdDetection();

        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_DIR));
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
                String sampleChordDirectory = mChordDir.replaceAll("\\*", sampleId);
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

            String samplePurpleDirectory = mPurpleDataDir.replaceAll("\\*", sampleId);
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
            String fileName = outputDir + "purple_hrd_analysis.csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,HRDStatus,HRD,BRCA1,BRCA2");
            writer.write(",LohSegments,SegmentBreaks,SegmentImbalances,Score");
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
            mWriter.write(format("%s,%s,%.3f,%.3f,%.3f",
                    sampleId, chordData.hrStatus(), chordData.hrdValue(), chordData.BRCA1Value(), chordData.BRCA2Value()));

            mWriter.write(format(",%d,%.1f,%d,%.1f",
                    hrdData.LohSegments, hrdData.SegmentBreaks, hrdData.SegmentImbalances, hrdData.score()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene overlap file output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addSampleIdFile(options);

        options.addOption(PURPLE_DATA_DIR, true, "Directory pattern for sample purple directory");
        addLoggingOptions(options);
        addOutputOptions(options);
        addThreadOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        HrdDetectionAnalyser hrdDetectionAnalyser = new HrdDetectionAnalyser(cmd);
        hrdDetectionAnalyser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }


}
