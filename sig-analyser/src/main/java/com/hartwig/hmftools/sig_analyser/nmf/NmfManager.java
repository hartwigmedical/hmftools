package com.hartwig.hmftools.sig_analyser.nmf;

import static java.lang.Integer.max;

import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.getNewFile;
import static com.hartwig.hmftools.common.sigs.SigMatrix.extractNonZeros;
import static com.hartwig.hmftools.common.sigs.SigMatrix.writeMatrixData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.SigMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfManager {

    private static final Logger LOGGER = LogManager.getLogger(NmfManager.class);

    private String mOutputDir;
    private String mOutputFileId;

    private GenericDataCollection mDataCollection;

    private SigMatrix mSampleCountsMatrix;
    private SigMatrix mReferenceSigs;
    private SigMatrix mReferenceContribs;

    private NmfCalculator mNmfCalculator;

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
                mNmfCalculator.setSignatures(mReferenceSigs);
        }

        if(!mConfig.RefContribFilename.isEmpty())
        {
            GenericDataCollection dataCollection = GenericDataLoader.loadFile(mConfig.RefContribFilename);
            mReferenceContribs = DataUtils.createMatrixFromListData(dataCollection.getData());
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

        SigMatrix contributions = sampleFitter.getContributions();
        contributions.cacheTranspose();

        writeSignatures(mReferenceSigs);
        writeContributions(contributions);

        mPerfCounter.logStats();
    }

    public void writeSignatures(final SigMatrix signatures)
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

    public void writeContributions(final SigMatrix contributions)
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

}
