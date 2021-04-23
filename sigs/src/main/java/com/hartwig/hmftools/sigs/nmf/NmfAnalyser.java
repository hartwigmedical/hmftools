package com.hartwig.hmftools.sigs.nmf;

import static java.lang.Integer.max;

import static com.hartwig.hmftools.common.utils.MatrixUtils.createMatrixFromListData;
import static com.hartwig.hmftools.common.utils.MatrixUtils.writeMatrixData;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sigs.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getNewFile;
import static com.hartwig.hmftools.common.utils.Matrix.extractNonZeros;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.Matrix;

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

public class NmfAnalyser
{
   private static final Logger LOGGER = LogManager.getLogger(NmfAnalyser.class);

    private String mOutputDir;
    private String mOutputFileId;

    private GenericDataCollection mDataCollection;

    private Matrix mSampleCountsMatrix;
    private Matrix mReferenceSigs;
    private Matrix mReferenceContribs;

    private NmfCalculator mNmfCalculator;

    private List<NmfRun> mRuns;

    private NmfConfig mConfig;

    PerformanceCounter mPerfCounter;

    public NmfAnalyser()
    {
        mOutputDir = "";
        mOutputFileId = "";
        mDataCollection = null;
        mSampleCountsMatrix = null;
        mReferenceSigs = null;
        mNmfCalculator = null;

        mRuns = Lists.newArrayList();

        mConfig = null;

        mPerfCounter = new PerformanceCounter("NMF");
    }

    public void initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mConfig = new NmfConfig(cmd);
        mDataCollection = collection;

        mPerfCounter.start("DataLoad");

        mSampleCountsMatrix = createMatrixFromListData(mDataCollection.getData());
        mSampleCountsMatrix = extractNonZeros(mSampleCountsMatrix);
        mSampleCountsMatrix.cacheTranspose();

        mNmfCalculator = new NmfCalculator(mSampleCountsMatrix, mConfig);

        if(!mConfig.RefSigFilename.isEmpty())
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(mConfig.RefSigFilename);
            mReferenceSigs = createMatrixFromListData(dataCollection.getData());
            mReferenceSigs.cacheTranspose();

            if(mConfig.UseRefSigs)
                mNmfCalculator.setSignatures(mReferenceSigs);
        }

        if(!mConfig.RefContribFilename.isEmpty())
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(mConfig.RefContribFilename);
            mReferenceContribs = createMatrixFromListData(dataCollection.getData());
            mReferenceContribs.cacheTranspose();

            if(!mConfig.FitRestrictToContribs)
            {
                mNmfCalculator.setContributions(mReferenceContribs);
            }
        }

        mPerfCounter.stop();
    }

    public void run() {

        if(mConfig.FitOnly)
            runFitOnly();
        else
            runNmf();
    }

    private void runNmf()
    {
        mPerfCounter.start("NMF");

        int startSigCount = mConfig.SigCount;
        int maxSigCount = max(mConfig.SigExpansionCount, mConfig.SigCount);

        double lowestRunScore = -1;
        int lowestRunIndex = -1;

        for (int sigCount = startSigCount; sigCount <= maxSigCount; ++sigCount)
        {
            LOGGER.info("starting run with sigCount({})", sigCount);

            NmfRun nmfRun = new NmfRun(mConfig, sigCount, mNmfCalculator, mReferenceSigs);

            if (!nmfRun.run()) {
                LOGGER.warn("run with sigCount({}) invalid, exiting", sigCount);
                break;
            }

            if (lowestRunScore < 0 || nmfRun.getLowestRunScore() < lowestRunScore)
            {
                lowestRunScore = nmfRun.getLowestRunScore();
                lowestRunIndex = mRuns.size();
            }

            mRuns.add(nmfRun);
        }

        mPerfCounter.stop();

        if (!mRuns.isEmpty())
        {
            final NmfRun nmfRun = mRuns.get(lowestRunIndex);

            if (mRuns.size() > 1)
            {
                LOGGER.info("optimal sigCount({})", nmfRun.getSigCount());
            }

            writeSignatures(nmfRun.getBestSignatures());
            writeContributions(nmfRun.getBestContributions());
        }

        mPerfCounter.logStats();
    }

    private void runFitOnly()
    {
        if(mReferenceSigs == null)
        {
            // for now only works with external sigs, in time could work with sig finder's set
            return;
        }

        mPerfCounter.start("NMF");

        NmfSampleFitter sampleFitter = new NmfSampleFitter(mConfig, mSampleCountsMatrix, mReferenceSigs);

        if(mConfig.FitRestrictToContribs)
        {
            sampleFitter.setRefContributions(mReferenceContribs);
        }

        sampleFitter.fitSamples();

        mPerfCounter.stop();

        if(!sampleFitter.isValid())
            return;

        Matrix contributions = sampleFitter.getContributions();
        contributions.cacheTranspose();

        writeSignatures(mReferenceSigs);
        writeContributions(contributions);

        mPerfCounter.logStats();
    }

    public void writeSignatures(final Matrix signatures)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir,mOutputFileId + "_nmf_sigs.csv");

            int i = 0;
            for(; i < signatures.Cols-1; ++i)
            {
                writer.write(String.format("%d,", i));
            }
            writer.write(String.format("%d", i));

            writer.newLine();

            writeMatrixData(writer, signatures, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeContributions(final Matrix contributions)
    {
        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_nmf_contribs.csv");

            final List<String> sampleNames = mDataCollection.getFieldNames();

            int i = 0;
            for(; i < sampleNames.size()-1; ++i)
            {
                writer.write(String.format("%s,", sampleNames.get(i)));
            }
            writer.write(String.format("%s", sampleNames.get(i)));

            writer.newLine();

            writeMatrixData(writer, contributions, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(SAMPLE_COUNTS_FILE, true, "Path to the main input file");
        options.addOption(OUTPUT_DIR, true, "Path to output files");

        NmfConfig.addCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        SIG_LOGGER.info("running NMF signature analyser");

        final GenericDataCollection collection = GenericDataLoader.loadFile(cmd.getOptionValue(SAMPLE_COUNTS_FILE));

        NmfAnalyser nmfAnalyser = new NmfAnalyser();
        nmfAnalyser.initialise(collection, cmd);
        nmfAnalyser.run();

        SIG_LOGGER.info("NMF analysis complete");
    }

}
