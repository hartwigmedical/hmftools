package com.hartwig.hmftools.linx.cn;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvaConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.linx.types.SvaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.types.SvaConfig.formOutputPath;

import java.io.BufferedWriter;
import java.io.IOException;

import java.util.List;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CopyNumberAnalyser
{

    private boolean mWriteAdjustedPloidyToFile;
    private boolean mWriteVerbosePloidyData;
    private boolean mWriteLohData;

    private List<String> mSampleIds;
    private BufferedWriter mFileWriter;
    private BufferedWriter mRecalcPloidyFileWriter;

    private static int PLOIDY_CALC_COLUMN_COUNT = 4;

    private static final String WRITE_PLOIDY_TO_FILE = "write_ploidy_data";
    private static final String WRITE_VERBOSE_PLOIDY_DATA = "verbose_ploidy_data";

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyser.class);

    public CopyNumberAnalyser(final String purpleDataPath, final String outputPath, DatabaseAccess dbAccess)
    {
        mFileWriter = null;
        mRecalcPloidyFileWriter = null;

        mWriteLohData = false;
        mWriteAdjustedPloidyToFile = false;
        mWriteVerbosePloidyData = false;

    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(WRITE_PLOIDY_TO_FILE, false, "Write adjusted ploidy to CSV");
        options.addOption(WRITE_VERBOSE_PLOIDY_DATA, false, "Write all ploidy calc working data");
    }

    public boolean loadConfig(final CommandLine cmd, final List<String> sampleIds)
    {
        mSampleIds.addAll(sampleIds);

        mWriteLohData = true;
        mWriteAdjustedPloidyToFile = cmd.hasOption(WRITE_PLOIDY_TO_FILE);
        mWriteVerbosePloidyData = cmd.hasOption(WRITE_VERBOSE_PLOIDY_DATA);

        return true;
    }


    private void loadSVData(final String sampleId)
    {
        /*
        mSvDataList = mDbAccess.readStructuralVariantData(sampleId);

        if(mSvDataList.isEmpty())
        {
            LOGGER.warn("sample({}) no SV records found", sampleId);
        }
        */
    }

    public void runAnalysis()
    {
        /*

                if (mWriteAdjustedPloidyToFile)
        {
            initialisePloidyWriter();
        }


        // mode for CN analyis in isolation from SV analysis
        if(mDbAccess == null)
        {
            LOGGER.warn("batch mode requires DB connection");
            return;
        }

        int sampleCount = 0;

        for(final String sampleId : mSampleIds)
        {
            LOGGER.info("analysing sample({}), totalProcessed({})", sampleId, sampleCount);

            loadSVData(sampleId);

            loadCopyNumberData(sampleId);

            processSampleData(sampleId);
            ++sampleCount;
        }

        */
    }

    /*
    private void processSampleData(final String sampleId)
    {
        linkCopyNumberAndSvData(sampleId);

        findLohEvents(sampleId);

        if(mWriteAdjustedPloidyToFile)
            reaclcAdjustedPloidy(sampleId);
    }
    */

    /*

            try
        {
            BufferedWriter writer = mRecalcPloidyFileWriter;

                if (writer != null)
                {
                    if (!mWriteVerbosePloidyData)
                    {
                        writer.write(String.format("%s,%d,%.4f,%.4f",
                                sampleId, svData.id(),
                                ploidyEstimate, ploidyUncertainty));
                    }
                    else
                    {
                        writer.write(String.format("%s,%d,%s,%.4f,%.4f,%.4f,%d,%d",
                                sampleId, svData.id(), svData.type(), svData.ploidy(), adjVafStart, adjVafEnd,
                                tumorReadCountStart, tumorReadCountEnd));

                        writer.write(String.format(",%s,%d,%d,%.4f,%.4f,%d,%d",
                                svData.startChromosome(), svData.startPosition(), svData.startOrientation(),
                                maxCNStart, cnChgStart, startDepthData[0], startDepthData[1]));

                        writer.write(String.format(",%s,%d,%d,%.4f,%.4f,%d,%d",
                                svData.endChromosome(), svData.endPosition(), svData.endOrientation(), maxCNEnd, cnChgEnd,
                                endDepthData != null ? endDepthData[0] : 0, endDepthData != null ? endDepthData[1] : 0));

                        writer.write(String.format(",%.2f,%.2f,%.2f,%.2f",
                                ploidyEstimate, ploidyUncertainty,
                                ploidyEstimate - ploidyUncertainty,
                                ploidyEstimate + ploidyUncertainty));
                    }

                    writer.newLine();
                }


     */


    private void writeLohData()
    {
        /*
                try
        {
            if(mWriteLohData)
            {
                if (mFileWriter == null)
                {
                    String outputFileName = mOutputPath;

                    outputFileName += "CN_LOH_EVENTS.csv";

                    mFileWriter = createBufferedWriter(outputFileName, false);

                    // SV info
                    mFileWriter.write("SampleId,Chromosome,CnIdStart,CnIdEnd,PosStart,PosEnd,SegStart,SegEnd,");
                    mFileWriter.write("PrevCN,StartCN,EndCN,MinCN,SegCount,Length,StartSV,EndSV,Skipped,IsValid");
                    mFileWriter.newLine();
                }
            }


                   if(mWriteLohData && mFileWriter != null)
        {
            mFileWriter.write(String.format("%s,%s,%d,%d,%d,%d,%s,%s",
                    sampleId, lohData.Chromosome, lohData.CnIdStart, lohData.CnIdEnd,
                    lohData.PosStart, lohData.PosEnd, lohData.SegStart, lohData.SegEnd));

            mFileWriter.write(String.format(",%.4f,%.4f,%.4f,%.4f,%d,%d,%s,%s,%s,%s",
                    lohData.PrevCN, lohData.StartCN, lohData.EndCN, lohData.MinCN,
                    lohData.SegCount, lohData.Length, lohData.StartSV, lohData.EndSV, lohData.Skipped, lohData.IsValid));

            mFileWriter.newLine();
        }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to copy number LOH outputFile: {}", e.toString());
            return 0;
        }

         */
    }

    private void initialisePloidyWriter()
    {
        try
        {
            if (mRecalcPloidyFileWriter == null)
            {
                String outputFileName = ""; // mOutputPath;

                outputFileName += "CN_PLOIDY_CALC_DATA.csv";

                mRecalcPloidyFileWriter = createBufferedWriter(outputFileName, false);

                if(!mWriteVerbosePloidyData)
                {
                    // SV info
                    mRecalcPloidyFileWriter.write("SampleId,SvId,EstPloidy,EstUncertainty");
                }
                else
                {
                    mRecalcPloidyFileWriter.write("SampleId,SvId,Type,Ploidy,VafStart,VafEnd,TumorRCStart,TumorRCEnd");
                    mRecalcPloidyFileWriter.write(",ChrStart,PosStart,OrientStart,MaxCNStart,CNChgStart,PrevDWCountStart,NextDWCountStart");
                    mRecalcPloidyFileWriter.write(",ChrEnd,PosEnd,OrientEnd,MaxCNEnd,CNChgEnd,PrevDWCountEnd,NextDWCountEnd");
                    mRecalcPloidyFileWriter.write(",EstPloidy,EstUncertainty,MinPloidy,MaxPloidy");
                }

                mRecalcPloidyFileWriter.newLine();
            }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to ploidy recalc outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
        closeBufferedWriter(mRecalcPloidyFileWriter);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String outputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));

        CopyNumberAnalyser cnAnalyser = new CopyNumberAnalyser("", outputDir, null);
        // cnAnalyser.loadConfig(cmd, samplesList);

        // run CN analysis, which will write a bunch of cohort-wide sample data, then exit
        cnAnalyser.runAnalysis();
        cnAnalyser.close();

        LOGGER.info("CN analysis complete");
    }


}
