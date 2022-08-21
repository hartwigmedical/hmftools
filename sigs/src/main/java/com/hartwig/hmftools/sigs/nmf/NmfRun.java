package com.hartwig.hmftools.sigs.nmf;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I1;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_I2;
import static com.hartwig.hmftools.sigs.common.CssRoutines.CSSR_VAL;
import static com.hartwig.hmftools.sigs.common.CssRoutines.getTopCssPairs;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sigs.common.SigReporter;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// A run is defined as a set of trials of differing starting conditions for a specific signature count

public class NmfRun {

    // record data for the best run
    private double mLowestResidualCount;
    private Matrix mBestSignatures;
    private Matrix mBestContributions;
    private List<Matrix> mUniqueSignatures;
    private Matrix mReferenceSigs; // for post-run sig comparison reporting
    private NmfCalculator mCalculator;

    private NmfConfig mConfig;

    private int mSigCount;
    private int mBucketCount;

    private boolean mValid;
    private Matrix mRandomStartSignatures;

    PerformanceCounter mPerfCounter;

    private static final Logger LOGGER = LogManager.getLogger(NmfRun.class);

    public NmfRun(final NmfConfig config, int sigCount, NmfCalculator nmfCalculator, final Matrix referenceSigs)
    {
        mConfig = config;
        mSigCount = sigCount;

        mCalculator = nmfCalculator;
        mValid = false;

        mBucketCount = mCalculator.getSampleCounts().Rows;

        mLowestResidualCount = -1;
        mBestSignatures = null;
        mBestContributions = null;
        mUniqueSignatures = Lists.newArrayList();

        mReferenceSigs = referenceSigs;

        mRandomStartSignatures = null;
        generateRandomStartSignatures();

        mPerfCounter = new PerformanceCounter(String.format("NMF %d Sigs", mSigCount));
    }

    public int getSigCount() { return mSigCount; }

    public boolean isValid() { return mValid; }

    public double getLowestRunScore() { return mLowestResidualCount; }

    public final Matrix getBestSignatures() { return mBestSignatures; }

    public final Matrix getBestContributions() { return mBestContributions; }

    public boolean run()
    {
        mValid = true;

        mCalculator.setSigCount(mSigCount);
        mCalculator.setRandomSignatures(mRandomStartSignatures);
        mCalculator.clearLowestCost();

        mPerfCounter.start();

        PerformanceCounter runPC = new PerformanceCounter("NMF Runs");

        boolean hasValidRun = false;

        for (int i = 0; i < mConfig.RunCount; ++i)
        {
            runPC.start();
            mCalculator.performRun(i);
            runPC.stop();

            if (!mCalculator.isValid()) {
                LOGGER.warn("exiting on invalid NMF run");
                mValid = false;
                break;
            }

            double newRunScore = mCalculator.getTotalResiduals();
            final Matrix newSigs = mCalculator.getSignatures();

            if (i == 0 || !hasValidRun)
            {
                hasValidRun = true;

                mLowestResidualCount = newRunScore;
                mBestSignatures = new Matrix(newSigs);
                mBestContributions = new Matrix(mCalculator.getContributions());

                // mUniqueSignatures.add(new SigMatrix(newSigs));
            }
            else
            {
                if (newRunScore < mLowestResidualCount)
                {
                    double residualsPercent = newRunScore / mCalculator.getTotalCount();
                    LOGGER.debug(String.format("run %d: score lowered(%.0f > %.0f) percent(%.5f)",
                            i, mLowestResidualCount, newRunScore, residualsPercent));

                    mLowestResidualCount = newRunScore;
                    mBestSignatures.setData(newSigs.getData());
                    mBestContributions.setData(mCalculator.getContributions().getData());
                }

                // store if this new signature is significantly different
                // cacheUniqueSignatures(newSigs); // currently unused
            }
        }

        mPerfCounter.stop();

        if (!mValid)
            return false;

        double bestFitPercent = mLowestResidualCount / mCalculator.getTotalCount();

        LOGGER.info(String.format("sigCount(%d) %d run(s) complete, lowestResiduals(%.0f perc=%.5f)",
                mSigCount, mConfig.RunCount, mLowestResidualCount, bestFitPercent));

        mBestSignatures.cacheTranspose();
        mBestContributions.cacheTranspose();

        SigReporter sigReporter = new SigReporter(mCalculator.getSampleCounts(), mBestSignatures, mBestContributions,
                mCalculator.getRefSignatures(), mReferenceSigs, mConfig);

        sigReporter.runAnalysis();

        mPerfCounter.logStats();
        runPC.logStats();

        return mValid;
    }

    private void cacheUniqueSignatures(final Matrix newSigs) {
        if (mUniqueSignatures.size() >= 10)
            return;

        boolean matchFound = false;
        for (final Matrix sig : mUniqueSignatures) {
            if (sig.equals(newSigs))
                continue;

            if (signaturesEqual(sig, newSigs)) {
                matchFound = true;
                break;
            }
        }

        if (!matchFound) {
            LOGGER.debug(String.format("storing new unique signature"));
            mUniqueSignatures.add(new Matrix(newSigs));
        }
    }

    public static boolean signaturesEqual(final Matrix sigs1, final Matrix sigs2) {
        // use CSS to compare each pair of sigs from the 2 sets
        // return true if the set of sigs are a close match
        double cssMatchCutoff = 0.98;
        List<double[]> cssResults = getTopCssPairs(sigs1, sigs2, cssMatchCutoff, true, false);

        int sigCount = sigs1.Cols;

        if (cssResults.size() != sigCount) {
            // LOGGER.debug("matched CSS sigCount({}) less than sigCount({})", cssResults.size(), sigCount);
            return false;
        }

        return true;
    }

    private void generateRandomStartSignatures() {
        // create a set of random bucket ratios for use in signatures
        // based on the frequency of bucket counts
        Random random = new Random(123456);

        int randomSigCount = 100; // could use combination of run count and sig count
        mRandomStartSignatures = new Matrix(mBucketCount, randomSigCount);

        // create a poisson distribution around the bucket ratio for the cohort
        final double[] bucketTotals = mCalculator.getBucketTotals();
        double totalCount = mCalculator.getTotalCount();

        double[][] rData = mRandomStartSignatures.getData();

        // double varianceFactor = 10;

        // assign bucket ratios randomly but favour buckets with higher cohort percentages
        for (int i = 0; i < bucketTotals.length; ++i) {

            double bucketPerc = bucketTotals[i] / totalCount;

            for (int j = 0; j < randomSigCount; ++j) {

                // double randomRatio = exp(log(bucketPerc) + (random.nextDouble() - 0.5) * varianceFactor);
                double randomRatio = sqrt(bucketPerc) * random.nextDouble();
                rData[i][j] = min(randomRatio, 1.0);
            }
        }

        // convert these to signatures (defined as ratios adding to 1)
        // and remove any that aren't sufficiently unique
        int sigIndex = 0;

        for (int i = 0; i < randomSigCount; ++i) {

            double bucketRatioTotal = 0;

            // for each signature, either set a random bucket ratio or take the ref one
            // and add up all ratios, so they can be used to split the actual bucket count amongst them
            for (int j = 0; j < mBucketCount; ++j) {

                bucketRatioTotal += rData[j][i];
            }

            for (int j = 0; j < mBucketCount; ++j) {

                rData[j][i] /= bucketRatioTotal;
            }
        }

        // check for similarities
        List<double[]> similarSigs = getTopCssPairs(mRandomStartSignatures, mRandomStartSignatures, 0.98, false, true);

        for(final double[] result : similarSigs)
        {
            LOGGER.debug(String.format("similar random sigs(%.0f & %.0f) have css(%.6f)", result[CSSR_I1], result[CSSR_I2], result[CSSR_VAL]));
        }
    }

}

