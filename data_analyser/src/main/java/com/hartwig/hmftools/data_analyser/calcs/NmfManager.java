package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Integer.max;

import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;
import static com.hartwig.hmftools.data_analyser.types.NmfMatrix.extractNonZeros;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfManager {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private String mOutputDir;
    private String mOutputFileId;

    private GenericDataCollection mDataCollection;

    private NmfMatrix mSampleCountsMatrix;
    private NmfMatrix mReferenceSigs;

    private NmfCalculator mNmfCalculator;
    private SigFinder mSigFinder;

    private List<NmfRun> mRuns;

    private NmfConfig mConfig;

    PerformanceCounter mPerfCounter;

    public NmfManager()
    {
        mOutputDir = "";
        mOutputFileId = "";
        mDataCollection = null;
        mSampleCountsMatrix = null;
        mReferenceSigs = null;
        mNmfCalculator = null;
        mSigFinder = null;

        mRuns = Lists.newArrayList();

        mConfig = null;

        mPerfCounter = new PerformanceCounter("NMF");
    }

    public void initialise(GenericDataCollection collection, final CommandLine cmd)
    {
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mConfig = new NmfConfig(cmd);
        mDataCollection = collection;

        mPerfCounter.start("DataLoad");

        mSampleCountsMatrix = DataUtils.createMatrixFromListData(mDataCollection.getData());
        mSampleCountsMatrix = extractNonZeros(mSampleCountsMatrix);
        mSampleCountsMatrix.cacheTranspose();

        mNmfCalculator = new NmfCalculator(mSampleCountsMatrix, mConfig);

        if(!mConfig.RefSigFilename.isEmpty())
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(mConfig.RefSigFilename);
            mReferenceSigs = DataUtils.createMatrixFromListData(dataCollection.getData());
            mReferenceSigs.cacheTranspose();

            if(mConfig.UseRefSigs)
                mNmfCalculator.setSignatures(mReferenceSigs, mConfig.SigFloatRate);
        }

        if(!mConfig.RefContribFilename.isEmpty())
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(mConfig.RefContribFilename);
            mNmfCalculator.setContributions(DataUtils.createMatrixFromListData(dataCollection.getData()));
        }

        mPerfCounter.stop();

        if(mConfig.FindSignatures)
        {
            mSigFinder = new SigFinder(mSampleCountsMatrix, mConfig, mOutputDir, mOutputFileId, mDataCollection.getFieldNames());

            mPerfCounter.start("SigFinder");
            mSigFinder.findSignatures();
            mPerfCounter.stop();

            if(mSigFinder.getSignatures() != null)
            {
                mNmfCalculator.setSignatures(mSigFinder.getSignatures(), mConfig.SigFloatRate);
            }

            if(mSigFinder.getContributions() != null)
            {
                mNmfCalculator.setContributions(mSigFinder.getContributions());
            }
        }
    }

    public void run() {
        mPerfCounter.start("NMF");

        int startSigCount = mConfig.SigCount;
        int maxSigCount = max(mConfig.SigExpansionCount, mConfig.SigCount);

        double lowestRunScore = -1;
        int lowestRunIndex = -1;

        for (int sigCount = startSigCount; sigCount <= maxSigCount; ++sigCount) {
            LOGGER.info("starting run with sigCount({})", sigCount);

            NmfRun nmfRun = new NmfRun(mConfig, sigCount, mNmfCalculator, mReferenceSigs);

            if (!nmfRun.run()) {
                LOGGER.warn("run with sigCount({}) invalid, exiting", sigCount);
                break;
            }

            if (lowestRunScore < 0 || nmfRun.getLowestRunScore() < lowestRunScore) {
                lowestRunScore = nmfRun.getLowestRunScore();
                lowestRunIndex = mRuns.size();
            }

            mRuns.add(nmfRun);
        }

        mPerfCounter.stop();

        if (!mRuns.isEmpty()) {

            final NmfRun nmfRun = mRuns.get(lowestRunIndex);

            if (mRuns.size() > 1) {
                LOGGER.info("optimal sigCount({})", nmfRun.getSigCount());
            }

            writeSignatures(nmfRun.getBestSignatures());
            writeContributions(nmfRun.getBestContributions());
        }

        mPerfCounter.logStats();

    }

    public void writeSignatures(final NmfMatrix signatures)
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
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeContributions(final NmfMatrix contributions)
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
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void runTests() {

        testMatrixCopy();
    }

    private void testMatrixCopy()
    {
        NmfMatrix m1 = new NmfMatrix(5, 5);

        Random random = new Random(123456);
        DataUtils.initRandom(m1, 0, 1, random);

        NmfMatrix m2 = new NmfMatrix(m1);

        m1.set(0, 0, 100);
    }

    private void testSigCompare()
    {
        NmfMatrix sigs1 = new NmfMatrix(5, 5);
        NmfMatrix sigs2 = new NmfMatrix(5, 5);

        double[][] s1Data = sigs1.getData();
        double[][] s2Data = sigs2.getData();

        for(int i = 0; i < sigs1.Rows; ++i)
        {
            for(int j = 0; j < sigs1.Cols; ++j)
            {
                s1Data[i][j] = (i+1) * (j+1);

                if(j < sigs1.Cols-1)
                    s2Data[i][j] = s1Data[i][j];
            }
        }

        NmfRun.signaturesEqual(sigs1, sigs2);
    }

}
