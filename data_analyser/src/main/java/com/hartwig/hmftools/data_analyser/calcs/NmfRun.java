package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I1;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_I2;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.CSSR_VAL;
import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.getTopCssPairs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// A run is defined as a set of trials of differing starting conditions for a specific signature count

public class NmfRun {

    // record data for the best run
    private double mLowestResidualCount;
    private NmfMatrix mBestSignatures;
    private NmfMatrix mBestContributions;
    private double[] mLowestSampleResiduals;
    private double[] mLowestBucketResiduals;
    private List<NmfMatrix> mUniqueSignatures;

    private NmfCalculator mCalculator;

    private NmfConfig mConfig;

    private int mSigCount;
    private int mSampleCount;
    private int mBucketCount;

    private boolean mValid;

    PerformanceCounter mPerfCounter;

    private static final Logger LOGGER = LogManager.getLogger(NmfRun.class);

    public NmfRun(final NmfConfig config, int sigCount, NmfCalculator nmfCalculator)
    {
        mConfig = config;
        mSigCount = sigCount;

        mCalculator = nmfCalculator;
        mValid = false;

        mSampleCount = mCalculator.getSampleCounts().Cols;
        mBucketCount = mCalculator.getSampleCounts().Rows;
        mLowestSampleResiduals = new double[mSampleCount];
        mLowestBucketResiduals = new double[mBucketCount];

        mLowestResidualCount = -1;
        mBestSignatures = null;
        mBestContributions = null;
        mUniqueSignatures = Lists.newArrayList();

        mPerfCounter = new PerformanceCounter(String.format("NMF %d Sigs", mSigCount));
    }

    public int getSigCount() { return mSigCount; }
    public boolean isValid() { return mValid; }
    public double getLowestRunScore() { return mLowestResidualCount; }
    public final NmfMatrix getBestSignatures() { return mBestSignatures; }
    public final NmfMatrix getBestContributions() { return mBestContributions; }

    public boolean run()
    {
        mValid = true;

        mCalculator.setSigCount(mSigCount);
        mCalculator.clearLowestCost();

        mPerfCounter.start();

        PerformanceCounter runPC = new PerformanceCounter("NMF Runs");

        boolean hasValidRun = false;

        for(int i = 0; i < mConfig.RunCount; ++i)
        {
            runPC.start();
            mCalculator.performRun();
            runPC.stop();

            if(!mCalculator.isValid())
            {
                LOGGER.warn("exiting on invalid NMF run");
                mValid = false;
                break;
            }

            double newRunScore = mCalculator.getTotalResiduals();
            final NmfMatrix newSigs = mCalculator.getSignatures();

            if(i == 0 || !hasValidRun)
            {
                hasValidRun = true;

                mLowestResidualCount = newRunScore;
                mBestSignatures = new NmfMatrix(newSigs);
                mBestContributions = new NmfMatrix(mCalculator.getContributions());
                cacheResidualData();

                // mUniqueSignatures.add(new NmfMatrix(newSigs));
            }
            else
            {
                if(newRunScore < mLowestResidualCount)
                {
                    LOGGER.debug(String.format("run %d: score lowered(%.1f > %.1f)", i, mLowestResidualCount, newRunScore));

                    mLowestResidualCount = newRunScore;
                    mBestSignatures.setData(newSigs.getData());
                    mBestContributions.setData(mCalculator.getContributions().getData());
                    cacheResidualData();
                }

                // store if this new signature is significantly different
                // cacheUniqueSignatures(newSigs); // currently unused
            }
        }

        mPerfCounter.stop();

        if(!mValid)
            return false;

        double bestFitPercent = mLowestResidualCount / mCalculator.getTotalCount();

        LOGGER.info(String.format("sigCount(%d) %d run(s) complete, lowestResiduals(%.0f perc=%.3f)",
                mSigCount, mConfig.RunCount, mLowestResidualCount, bestFitPercent));

        logResidualData();
        logSignatureSummaryData();

        mPerfCounter.logStats();
        runPC.logStats(false); // mConfig.LogVerbose

        return mValid;
    }

    private void cacheResidualData()
    {
        copyVector(mCalculator.getSampleResiduals(), mLowestSampleResiduals);
        copyVector(mCalculator.getBucketResiduals(), mLowestBucketResiduals);
    }

    private void logResidualData()
    {
        // report any stand-out sample or bucket, defined as X times the average
        double residualFactor = 5;

        double totalCount = mCalculator.getTotalCount();
        double totalResidualsPerc = mLowestResidualCount / totalCount;
        double resPercThreshold = totalResidualsPerc * residualFactor;
        double sampleCountThreshold = mCalculator.getTotalCount() / (double)mSampleCount * 0.25;

        for(int i = 0; i < mSampleCount; ++i) {

            double resCount = mLowestSampleResiduals[i];
            double sampleCount = mCalculator.getSampleTotals()[i];

            if (resCount == 0)
                continue;

            double resPerc = resCount / sampleCount;

            if (resPerc >= resPercThreshold && sampleCount >= sampleCountThreshold) {

                LOGGER.debug(String.format("sample(%d) residuals(count=%.0f perc=%.2f) above avg(%.3f)",
                        i, resCount, resPerc, totalResidualsPerc));
            }
        }

        double bucketCountThreshold = totalCount / (double)mBucketCount * 0.25;

        for(int i = 0; i < mBucketCount; ++i) {

            double resCount = mLowestBucketResiduals[i];
            double bucketCount = mCalculator.getBucketTotals()[i];

            if (resCount == 0)
                continue;

            double resPerc = resCount / bucketCount;

            if (resPerc >= resPercThreshold && bucketCount >= bucketCountThreshold) {

                LOGGER.debug(String.format("bucket(%d) residuals(count=%.0f perc=%.2f) above avg(%.3f)",
                        i, resCount, resPerc, totalResidualsPerc));
            }
        }
    }

    private void cacheUniqueSignatures(final NmfMatrix newSigs)
    {
        if(mUniqueSignatures.size() >= 10)
        return;

        boolean matchFound = false;
        for (final NmfMatrix sig : mUniqueSignatures) {
            if (sig.equals(newSigs))
                continue;

            if (signaturesEqual(sig, newSigs)) {
                matchFound = true;
                break;
            }
        }

        if (!matchFound) {
            LOGGER.debug(String.format("storing new unique signature"));
            mUniqueSignatures.add(new NmfMatrix(newSigs));
        }
    }

    public void logSignatureSummaryData()
    {
        if(mBestSignatures == null)
            return;

        CosineSim.logSimilarites(mBestSignatures, 0.98, "sig");

        if(mConfig.LogVerbose && mCalculator.getRefSignatures() != null && mConfig.SigFloatRate > 0)
        {
            // compare the starting ref sigs to the adjusted ones
            double cssMatchCutoff = 0.90;
            List<double[]> cssResults = getTopCssPairs(mCalculator.getRefSignatures(), mBestSignatures, cssMatchCutoff, true, false);

            if(cssResults.isEmpty())
            {
                LOGGER.debug("no similar sigs between ref and discovered");
            }
            else
            {
                for(final double[] result : cssResults)
                {
                    LOGGER.debug(String.format("ref sig(%.0f) matches adjusted sig(%.0f) with css(%.4f)",
                            result[CSSR_I1], result[CSSR_I2], result[CSSR_VAL]));
                }
            }
        }
    }

    public static boolean signaturesEqual(final NmfMatrix sigs1, final NmfMatrix sigs2)
    {
        // use CSS to compare each pair of sigs from the 2 sets
        // return true if the set of sigs are a close match
        double cssMatchCutoff = 0.98;
        List<double[]> cssResults = getTopCssPairs(sigs1, sigs2, cssMatchCutoff, true, false);

        int sigCount = sigs1.Cols;

        if(cssResults.size() != sigCount)
        {
            // LOGGER.debug("matched CSS sigCount({}) less than sigCount({})", cssResults.size(), sigCount);
            return false;
        }

        return true;
    }

}
