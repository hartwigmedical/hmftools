package com.hartwig.hmftools.purple.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.containsFlag;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.region.ObservedRegion;

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
    private final int mThreads;

    private final HrdDetection mHrdDetection;
    private final DatabaseAccess mDbAccess;

    private final BufferedWriter mWriter;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String PURPLE_DATA_DIR = "purple_data_dir";
    private static final String THREADS = "threads";

    public HrdDetectionAnalyser(final CommandLine cmd)
    {
        mSampleIds = loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE));
        mPurpleDataDir = cmd.getOptionValue(PURPLE_DATA_DIR);
        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        mDbAccess = DatabaseAccess.createDatabaseAccess(cmd);
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

        PPL_LOGGER.info("Purple HRD analysis cnmplete");
    }

    private void processSample(final String sampleId)
    {
        final List<PurpleCopyNumber> copyNumbers = mDbAccess.readCopynumbers(sampleId);
        final PurityContext purityContext = mDbAccess.readPurityContext(sampleId);
        final ChordAnalysis chordAnalysis = mDbAccess.readChord(sampleId);

        if(chordAnalysis == null || purityContext == null)
        {
            PPL_LOGGER.info("sample({}) invalid purity({}) or chord({}) data",
                    sampleId, purityContext != null ? "valid" : "missing", chordAnalysis != null ? "valid" : "missing");
            return;
        }

        PPL_LOGGER.debug(format("sample(%s) ploidy(%.1f) cnRecords(%d) chord(%s %.3f)",
                sampleId, purityContext.bestFit().ploidy(), copyNumbers.size(), chordAnalysis.hrStatus(), chordAnalysis.hrdValue()));

        final HrdData hrdData = mHrdDetection.calculateHrdData(copyNumbers, purityContext.bestFit().ploidy());

        writeSampleData(sampleId, chordAnalysis, hrdData);
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

    private synchronized void writeSampleData(final String sampleId, final ChordAnalysis chordAnalysis, final HrdData hrdData)
    {
        try
        {
            mWriter.write(format("%s,%s,%.3f,%.3f,%.3f",
                    sampleId, chordAnalysis.hrStatus(), chordAnalysis.hrdValue(), chordAnalysis.BRCA1Value(), chordAnalysis.BRCA2Value()));

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
        options.addOption(SAMPLE_ID_FILE, true, "Sample ID file");
        options.addOption(THREADS, true, "Thread count, default = 0 (disabled)");

        DatabaseAccess.addDatabaseCmdLineArgs(options);
        options.addOption(PURPLE_DATA_DIR, true, "Directory pattern for sample purple directory");
        addLoggingOptions(options);
        addOutputOptions(options);

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
