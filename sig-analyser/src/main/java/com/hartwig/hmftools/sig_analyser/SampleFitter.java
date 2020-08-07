package com.hartwig.hmftools.sig_analyser;

import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigUtils.RESIDUAL_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigUtils.RESIDUAL_PERC;
import static com.hartwig.hmftools.common.sigs.SigUtils.RESIDUAL_TOTAL;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleListFile;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleMatrixCounts;
import static com.hartwig.hmftools.sig_analyser.loaders.SigDataLoader.DB_URL;
import static com.hartwig.hmftools.sig_analyser.loaders.SigDataLoader.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.sig_analyser.loaders.SigDataLoader.databaseAccess;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sig_analyser.loaders.SigSnvLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SampleFitter
{
    // config
    private static final String SIGNATURES_FILE = "mSignatures_file";
    private static final String SAMPLE_IDS = "sample";
    
    private final String mSampleIdsConfig;
    private final List<String> mSampleIdList;
    private final String mSnvCountsFile;
    private final String mSignaturesFile;
    private final String mOutputDir;

    private SigMatrix mSampleCountsMatrix;
    private SigMatrix mSignatures;
    private DatabaseAccess mDbAccess;
    private BufferedWriter mFitWriter;

    public SampleFitter(final CommandLine cmd)
    {
        mSnvCountsFile = cmd.getOptionValue(SAMPLE_COUNTS_FILE);
        mSignaturesFile = cmd.getOptionValue(SIGNATURES_FILE);
        mSampleIdsConfig = cmd.getOptionValue(SAMPLE_IDS);
        mSampleIdList = Lists.newArrayList();

        mSampleCountsMatrix = null;
        mSignatures = null;

        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if(!outputDir.endsWith(File.separator))
        {
            outputDir += File.separator;
        }

        mOutputDir = outputDir;
        mFitWriter = null;
        mDbAccess = null;

        if(cmd.hasOption(DB_URL))
        {
            try
            {
                mDbAccess = databaseAccess(cmd);

            } catch (SQLException e)
            {
                SIG_LOGGER.error("DB connection failed: {}", e.toString());
            }
        }
    }

    public void run()
    {
        SIG_LOGGER.info("running sample signature fit");

        if(!initialise())
            return;

        SIG_LOGGER.info("fitting {} samples");

        performFit();

        SIG_LOGGER.info("sample signature fit complete");

        closeBufferedWriter(mFitWriter);
    }
    
    private boolean initialise()
    {
        if(mSignaturesFile == null)
        {
            SIG_LOGGER.error("missing mSignatures file");
            return false;
        }

        if(mSampleIdsConfig == null && mSnvCountsFile == null)
        {
            SIG_LOGGER.error("missing sampleIds config and sample count file");
            return false;
        }

        final GenericDataCollection sigsCollection = GenericDataLoader.loadFile(mSignaturesFile);
        mSignatures = DataUtils.createMatrixFromListData(sigsCollection.getData());

        if(mSnvCountsFile != null)
        {
            mSampleCountsMatrix = loadSampleMatrixCounts(mSnvCountsFile, mSampleIdList);
        }
        else
        {
            // load from file or delimitered list
            if(mSampleIdsConfig.contains(".csv"))
            {
                // load from file
                mSampleIdList.addAll(loadSampleListFile(mSampleIdsConfig));
            }
            else
            {
                mSampleIdList.add(mSampleIdsConfig);
            }

            loadSnvCountsData();
        }

        initialiseOutputFiles();
        return true;
    }

    private void loadSnvCountsData()
    {
        if(mDbAccess == null)
            return;

        SigSnvLoader snvLoader = new SigSnvLoader(null, mSampleIdList, mOutputDir);
        snvLoader.loadData(mDbAccess);

        final String filename = mSampleIdList.size() == 1 ? mSampleIdList.get(0) + ".sig.snv_counts.csv" : "SIG_SNV_COUNTS.csv";
        snvLoader.writeSampleCounts(filename);

        mSampleCountsMatrix = snvLoader.getSampleBucketCounts();
    }
    
    private void performFit()
    {
        int sampleCount = mSampleCountsMatrix.Cols;
        int sigCount = mSignatures.Cols;
        
        final SigMatrix sampleContribs = new SigMatrix(sigCount, sampleCount);

        SIG_LOGGER.info("fitting sample({}) with {} mSignatures", sampleCount, sigCount);

        LeastSquaresFit lsqFit = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);

        for(int i = 0; i < sampleCount; ++i)
        {
            final double[] sampleCounts = mSampleCountsMatrix.getCol(i);
            lsqFit.initialise(mSignatures.getData(), sampleCounts);
            lsqFit.solve();

            final double[] sigAllocs = lsqFit.getContribs();
            sampleContribs.setCol(i, sigAllocs);
        }

        // post-run analysis
        for(int i = 0; i < sampleCount; ++i)
        {
            final String sampleId = mSampleIdList.get(i);
            final double[] sampleCounts = mSampleCountsMatrix.getCol(i);
            final double[] sigAllocs = sampleContribs.getCol(i);

            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal == 0)
                continue;

            final double[] fittedCounts = calculateFittedCounts(mSignatures, sigAllocs);
            final double[] residuals = calcResiduals(sampleCounts, fittedCounts, sampleTotal);

            double allocTotal = sumVector(sigAllocs);
            double allocPerc = allocTotal / sampleTotal;

            SIG_LOGGER.debug(String.format("sample(%s) alloc(%s perc=%.3f of total=%s) residuals(%.3f total=%s excess=%s)",
                    sampleId, sizeToStr(allocTotal), allocPerc, sizeToStr(sampleTotal),
                    residuals[RESIDUAL_PERC], sizeToStr(residuals[RESIDUAL_TOTAL]), sizeToStr(residuals[RESIDUAL_EXCESS])));

            if(SIG_LOGGER.isDebugEnabled())
            {
                List<Integer> sortedSigs = getSortedVectorIndices(sigAllocs, false);
                for(Integer sigIndex : sortedSigs)
                {
                    double sigAlloc = sigAllocs[sigIndex];
                    double sigPercent = sigAlloc / sampleTotal;

                    SIG_LOGGER.trace(String.format("sample(%s) sampleTotal(%.0f) sig(%d) alloc(%.0f perc=%.3f)",
                            sampleId, sampleTotal, sigIndex, sigAlloc, sigPercent));

                    if(sigAlloc < 1)
                        break;
                }
            }
        }

        /*
        try
        {
            BufferedWriter writer = getNewFile(outputDir, outputFile);
            writeMatrixData(writer, scCollection.getFieldNames(), sampleContribs, false);
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write sample contributions: {}", e.toString());
            return;
        }

         */
    }

    private void initialiseOutputFiles()
    {
        final String filename = mSampleIdList.size() == 1 ?
                SignatureAllocationFile.generateFilename(mOutputDir, mSampleIdList.get(0)) : mOutputDir + "SIG_SNV_FIT.csv";

        try
        {
            mFitWriter = createBufferedWriter(filename, false);
            mFitWriter.write(SignatureAllocationFile.header());
            mFitWriter.newLine();
        }
        catch (final IOException e)
        {
            SIG_LOGGER.error("error writing to outputFile({}): {}", filename, e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(SAMPLE_COUNTS_FILE, true, "Path to the main input file");
        options.addOption(SIGNATURES_FILE, true, "Signature definitions");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        addDatabaseCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SampleFitter sampleFitter = new SampleFitter(cmd);

    }
}
