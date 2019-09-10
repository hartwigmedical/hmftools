package com.hartwig.hmftools.linx.cn;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.linx.LinxConfig.DB_PASS;
import static com.hartwig.hmftools.linx.LinxConfig.DB_URL;
import static com.hartwig.hmftools.linx.LinxConfig.DB_USER;
import static com.hartwig.hmftools.linx.LinxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.LinxConfig.databaseAccess;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;

import java.io.BufferedWriter;
import java.io.IOException;

import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.linx.types.SvVarData;
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

    private boolean mWriteVerbosePloidyData;

    private final String mOutputPath;
    private final DatabaseAccess mDbAccess;

    private final CnDataLoader mCnDataLoader;

    private List<String> mSampleIds;
    private BufferedWriter mFileWriter;
    private BufferedWriter mRecalcPloidyFileWriter;

    private static int PLOIDY_CALC_COLUMN_COUNT = 4;

    private static final String WRITE_PLOIDY_TO_FILE = "write_ploidy_data";
    private static final String WRITE_VERBOSE_PLOIDY_DATA = "verbose_ploidy_data";

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberAnalyser.class);

    public CopyNumberAnalyser(final String outputPath, DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mOutputPath = outputPath;

        mCnDataLoader = new CnDataLoader("", dbAccess);

        mFileWriter = null;
        mRecalcPloidyFileWriter = null;

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

        mWriteVerbosePloidyData = cmd.hasOption(WRITE_VERBOSE_PLOIDY_DATA);
        return true;
    }

    public void runAnalysis()
    {
        final List<String> samplesList = mDbAccess.getSampleIds();

        int sampleCount = 0;
        for (final String sampleId : samplesList)
        {
            List<StructuralVariantData> svRecords = mDbAccess.readStructuralVariantData(sampleId);

            if (svRecords.isEmpty())
            {
                continue;
            }

            LOGGER.info("analysing sample({}), totalProcessed({})", sampleId, sampleCount);

            mCnDataLoader.loadSampleData(sampleId, svRecords);

            writeLohData(sampleId);
            ++sampleCount;
        }
    }

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


    private void writeLohData(final String sampleId)
    {
        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputPath + "CN_LOH_EVENTS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Chromosome,PosStart,PosEnd,SegStart,SegEnd");
                mFileWriter.write(",SegCount,Length,StartSV,EndSV");
                mFileWriter.newLine();
            }

            final List<LohEvent> lohEvents = mCnDataLoader.getLohData();

            for(final LohEvent lohData : lohEvents)
            {
                mFileWriter.write(String.format("%s,%s,%d,%d,%s,%s",
                        sampleId, lohData.Chromosome, lohData.PosStart, lohData.PosEnd, lohData.SegStart, lohData.SegEnd));

                mFileWriter.write(String.format(",%d,%d,%s,%s",
                        lohData.SegCount, lohData.length(), lohData.StartSV, lohData.EndSV));

                mFileWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to copy number LOH outputFile: {}", e.toString());
        }
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
    }

    public static void main(@NotNull final String[] args) throws ParseException, SQLException
    {
        final Options options = new Options();
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String outputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));

        final DatabaseAccess dbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        CopyNumberAnalyser cnAnalyser = new CopyNumberAnalyser(outputDir, dbAccess);
        // cnAnalyser.loadConfig(cmd, samplesList);

        cnAnalyser.runAnalysis();
        cnAnalyser.close();

        LOGGER.info("CN analysis complete");
    }


}
