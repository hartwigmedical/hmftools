package com.hartwig.hmftools.sigs.common;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.common.CssRoutines.getTopCssPairs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.nmf.NmfConfig;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigReporter {

    final private int mBucketCount;
    final private int mSampleCount;
    final private int mSigCount;

    // primary input - bucket counts per sample
    final Matrix mSampleCounts;
    final Matrix mContributions;
    final Matrix mSignatures;
    final Matrix mStartSignatures; // from beginning of NMF run
    final Matrix mRefSignatures; // external sigs to compare with
    final Matrix mFittedCounts;
    final NmfConfig mConfig;
    final private double mTotalCount;
    final private double[] mBucketTotals;
    final private double[] mSampleTotals;

    private double mTotalResiduals;
    private double mNetResiduals;
    private double[] mSampleResiduals;
    private double[] mBucketResiduals;

    private static int TOP_BUCKET_COUNT = 5;
    private static double TOP_BUCKET_PERC = 0.9;

    private static final Logger LOGGER = LogManager.getLogger(SigReporter.class);

    public SigReporter(final Matrix sampleCounts, final Matrix signatures, final Matrix contributions,
            final Matrix startSigs, final Matrix refSigs, final NmfConfig config)
    {
        mSampleCounts = sampleCounts;
        mSignatures = signatures;
        mContributions = contributions;
        mStartSignatures = startSigs;
        mRefSignatures = refSigs;
        mConfig = config;

        mSampleCount = mContributions.Cols;
        mBucketCount = mSignatures.Rows;
        mTotalCount = mSampleCounts.sum();
        mSigCount = mSignatures.Cols;

        mBucketTotals = new double[mBucketCount];
        for(int i = 0; i < mSampleCounts.Rows; ++i)
        {
            mBucketTotals[i] = sumVector(mSampleCounts.getRow(i));
        }

        mSampleTotals = new double[mSampleCounts.Cols];

        for(int i = 0; i < mSampleCounts.Cols; ++i)
        {
            mSampleTotals[i] = sumVector(mSampleCounts.getCol(i));
        }

        mSampleResiduals = new double[mSampleCount];
        mBucketResiduals = new double[mBucketCount];
        mTotalResiduals = 0;
        mNetResiduals = 0;

        mFittedCounts = mSignatures.multiply(mContributions);
    }

    public void runAnalysis()
    {
        logSignatureData();
        compareSignatures();
        // compareContributions();
        calcResidualData();
        logResidualData();
    }

    private void logSignatureData()
    {
        int top5PercIndex = (int)round(mBucketCount * 0.05);
        int top20PercIndex = (int)round(mBucketCount * 0.2);
        double avgSigContrib = 1 / (double)mSigCount;

        for(int sigId = 0; sigId < mSigCount; ++sigId)
        {
            double[] sampleContribs = mContributions.getRow(sigId);
            double sigTotal = sumVector(sampleContribs);

            // count number of samples where this signature is above the avg signature contribution
            int countAboveAvg = 0;

            for(int j = 0; j < mSampleCount; ++j) {

                double sampleSigContrib = mContributions.get(sigId, j);
                double samplePerc = sampleSigContrib / mSampleTotals[j];

                if (samplePerc > avgSigContrib) {

                    //                    LOGGER.debug(String.format("sig(%d) sample(%d) sigPerc(%.3f) count(%.0f of %.0f)",
                    //                            sigId, j, samplePerc, sampleSigContrib, mSampleTotals[j]));
                    ++countAboveAvg;
                }
            }

            final List<Integer> sortedContribs = getSortedVectorIndices(sampleContribs, false);
            double maxContrib = sampleContribs[sortedContribs.get(0)];

            // take median from samples with non-zero contribs
            int positiveCount = 0;
            for(int i = 0; i < sortedContribs.size(); ++i)
            {
                if(sampleContribs[sortedContribs.get(i)] < 1)
                    break;

                ++positiveCount;
            }

            int medIndex = positiveCount / 2;
            double medContrib = sampleContribs[sortedContribs.get(medIndex)];

            LOGGER.debug(String.format("sig(%d) summary: total(%.0f perc=%.3f) stats(med=%.0f max=%.0f) sample over-allocated(%d perc=%.3f)",
                    sigId, sigTotal, sigTotal/mTotalCount, medContrib, maxContrib, countAboveAvg, countAboveAvg/(double)mSampleCount));

            // log top-X buckets
            final double[] bucketRatios = mSignatures.getCol(sigId);

            final List<Integer> sortedIndices = getSortedVectorIndices(bucketRatios, false);

            // create the fitted bucket counts for this signature
            double[] sigBucketContribs = new double[mBucketCount];

            for(int j = 0; j < mSampleCount; ++j) {

                for (int k = 0; k < mBucketCount; ++k) {

                    sigBucketContribs[k] += mContributions.get(sigId, j) * bucketRatios[k];
                }
            }

            double percTotal = 0;
            int bucketIndex = 0;
            double top5PercTotal = 0;
            double top20PercTotal = 0;

            LOGGER.debug("sig({}) top buckets", sigId);

            while(bucketIndex <= top20PercIndex)
            {
                int bucketId = sortedIndices.get(bucketIndex);
                double bucketPerc = bucketRatios[bucketId];

                if(bucketIndex < TOP_BUCKET_COUNT && percTotal < TOP_BUCKET_PERC) {

                    double allocatedCount = sigBucketContribs[bucketId];
                    double allocatedPerc = allocatedCount / mBucketTotals[bucketId];

                    LOGGER.debug(String.format("sig(%d) %d: bucket(%d) sigPerc(%.4f) allocated(%.0f ofBucket=%.3f)",
                            sigId, bucketIndex, bucketId, bucketPerc, allocatedCount, allocatedPerc));
                }

                if(bucketIndex <= top5PercIndex)
                    top5PercTotal += bucketPerc;

                top20PercTotal += bucketPerc;
                percTotal += bucketPerc;
                ++bucketIndex;
            }

            LOGGER.debug(String.format("sig(%d) cumulative bucket perc: top-5(%.3f) top-20(%.3f)", sigId, top5PercTotal, top20PercTotal));
        }
    }

    private void compareSignatures()
    {
        if(mConfig.FitOnly)
            return;

        // the CSS cutoff can be lower since only interested in any similarity
        double cssCutoff = 0.95;

        // first look for similar sigs out of the NMF process
        List<double[]> cssResults = CssRoutines.getTopCssPairs(mSignatures, mSignatures, cssCutoff, false, true);

        if(cssResults.isEmpty())
        {
            LOGGER.debug("no similarities between resultant signatures");
        }
        else
        {
            LOGGER.debug("{} similar signatures:", cssResults.size());

            for (final double[] result : cssResults)
            {
                int sig1 = (int)result[CssRoutines.CSSR_I1];
                int sig2 = (int)result[CssRoutines.CSSR_I2];

                LOGGER.debug(String.format("sigs(%d and %d) similar with css(%.4f)",
                        sig1, sig2, result[CssRoutines.CSSR_VAL]));

                // run another CSS check over all matching buckets
                int matchingBuckets = 0;

                final double[] buckets1 = mSignatures.getCol(sig1);
                final double[] buckets2 = mSignatures.getCol(sig2);

                final List<Integer> sortedIndices1 = getSortedVectorIndices(buckets1, false);
                final List<Integer> sortedIndices2 = getSortedVectorIndices(buckets2, false);

                int sIndex = 0;
                while(sortedIndices1.get(sIndex) == sortedIndices2.get(sIndex))
                {
                    ++matchingBuckets;
                    ++sIndex;
                }

                if(matchingBuckets < 2)
                    continue;

                double[] mb1 = new double[matchingBuckets];
                double[] mb2 = new double[matchingBuckets];

                for(int i = 0; i < matchingBuckets; ++i)
                {
                    mb1[i] = buckets1[sortedIndices1.get(i)];
                    mb2[i] = buckets2[sortedIndices1.get(i)];
                }

                double css = calcCosineSim(mb1, mb2);

                LOGGER.debug(String.format("sigs(%d and %d) matching buckets(%d) with css(%.4f)",
                        sig1, sig2, matchingBuckets, css));
            }
        }

        // compare the starting ref sigs to the adjusted ones
        if (mStartSignatures != null && mConfig.SigFloatRate > 0)
        {
            cssResults = CssRoutines.getTopCssPairs(mStartSignatures, mSignatures, 0.98, true, false);

            if (cssResults.isEmpty())
            {
                LOGGER.debug("no similar sigs between starting and adjusted");
            }
            else
            {
                for (final double[] result : cssResults) {
                    LOGGER.debug(String.format("starting sig(%.0f) matches adjusted sig(%.0f) with css(%.4f)",
                            result[CssRoutines.CSSR_I1], result[CssRoutines.CSSR_I2], result[CssRoutines.CSSR_VAL]));
                }
            }
        }

        // compare an externally provided set of sigs (eg from Cosmic) to the generated ones
        if (mRefSignatures != null && !mConfig.UseRefSigs)
        {
            cssResults = CssRoutines.getTopCssPairs(mRefSignatures, mSignatures, 0.8, true, false);

            if (cssResults.isEmpty())
            {
                LOGGER.debug("no similar sigs between external ref and adjusted");
            }
            else
            {
                for (final double[] result : cssResults) {
                    int externalSigId = (int)result[CssRoutines.CSSR_I1] + 1; // bumped up to correspond to convention of starting with 1
                    int denovoSigId = (int)result[CssRoutines.CSSR_I2] + 1;
                    LOGGER.debug(String.format("external ref sig(%d) matches adjusted sig(%d) with css(%.4f)",
                            externalSigId, denovoSigId, result[CssRoutines.CSSR_VAL]));
                }
            }
        }
    }

    private void compareContributions()
    {
        // look for samples which are assigned very similar ratios of signatures
        List<double[]> contribsResults = CssRoutines.getTopCssPairs(mContributions, mContributions, 0.995, false, true);

        if (contribsResults.isEmpty())
        {
            LOGGER.debug("no similar sample contributions");
        }
        else
        {
            double totalPairs = mSampleCount * (mSampleCount + 1) * 0.5;
            double similarPerc = contribsResults.size() / totalPairs;
            LOGGER.debug(String.format("%d similar sample contributions found (perc=%.2f):", contribsResults.size(), similarPerc));

            for (final double[] result : contribsResults)
            {
                LOGGER.debug(String.format("sample1(%.0f) sample2(%.0f) with css(%.4f)",
                        result[CssRoutines.CSSR_I1], result[CssRoutines.CSSR_I2], result[CssRoutines.CSSR_VAL]));
            }
        }

    }

    private static double RESIDUAL_PROB_CUTOFF = 0.001;

    // indices into the results array
    private static int SAMP_ID = 0;
    private static int BUCK_ID = 1;
    private static int PROB = 2;
    private static int SB_CONTRIB = 3;
    private static int BC = 4;

    private void calcResidualData()
    {
        final double[][] vData = mFittedCounts.getData();
        final double[][] scData = mSampleCounts.getData();

        // create a list of residual differences that are likely to be larger than noise would allow
        // index 1 = sample, index 2 = bucket, index 3 = poisson prob
        List<double[]> ppResults = Lists.newArrayList();
        int totalLowProbDiff = 0;
        int lowLLProbSampleCount = 0;

        for(int n = 0; n < mSampleCount; ++n)
        {
            mSampleResiduals[n] = 0;

            for(int b = 0; b < mBucketCount; ++b)
            {
                mBucketResiduals[b] = 0;
            }

            int bucketLowProbCount = 0;

            for(int b = 0; b < mBucketCount; ++b)
            {
                int bucketCount = (int)scData[b][n];
                int sbContrib = (int)round(vData[b][n]);

                double diff = bucketCount - sbContrib;
                double absDiff = abs(diff);

                mSampleResiduals[n] += absDiff;
                mBucketResiduals[b] += absDiff;
                mNetResiduals += diff;
                mTotalResiduals += absDiff;

                if(sbContrib <= 0)
                    continue;

                double pProb = 0;

                if(poissonProbEarlyExit(sbContrib, bucketCount))
                    continue;

                PoissonDistribution poisson = new PoissonDistribution(sbContrib);

                if(sbContrib > bucketCount)
                {
                    pProb = poisson.cumulativeProbability(bucketCount);
                }
                else
                {
                    pProb = 1 - poisson.cumulativeProbability(bucketCount - 1);
                }

                if(pProb >= RESIDUAL_PROB_CUTOFF)
                    continue;

                totalLowProbDiff += abs(sbContrib - bucketCount);
                ++bucketLowProbCount;

                // add in order of lowest probability first
                int index = 0;
                for(;index < ppResults.size(); ++index)
                {
                    final double[] result = ppResults.get(index);
                    if(pProb < result[PROB])
                        break;
                }

                double[] probResult = {n, b, pProb, sbContrib, bucketCount};
                ppResults.add(index, probResult);
            }

            // test differences against a log-likehood of poisson noise probability
            final double[] actualCounts = mSampleCounts.getCol(n);
            final double[] fittedCounts = mFittedCounts.getCol(n);

            double llProb = CssRoutines.calcLogLikelihood(actualCounts, fittedCounts, false);

            if(llProb < 0.95)
                ++lowLLProbSampleCount;

            double sampleResidualPerc = mSampleResiduals[n] / mSampleTotals[n];

            LOGGER.debug(String.format("sample(%d) residuals(%.0f vs act=%.0f perc=%.4f) bucketLowProbCount(%d) llProb(%.6f)",
                    n, mSampleResiduals[n], mSampleTotals[n], sampleResidualPerc, bucketLowProbCount, llProb));
        }

        LOGGER.debug(String.format("sample-bucket residuals with low-prob: count(%d) residuals(%d perc=%.3f)",
                ppResults.size(), totalLowProbDiff, totalLowProbDiff/mTotalResiduals));

//        for(int i = 0; i < min(ppResults.size(), 20); ++i)
//        {
//            final double[] result = ppResults.get(i);
//            LOGGER.debug(String.format("sample(%.0f) bucket(%.0f) prob(%.4g) counts(bc=%.0f fit=%.0f)",
//                    result[SAMP_ID], result[BUCK_ID], result[PROB], result[SB_CONTRIB], result[BC]));
//        }

        LOGGER.info(String.format("total residuals(%.0f vs count=%.0f perc=%.5f) low sample log-likelihood(%d perc=%.2f)",
                mTotalResiduals, mTotalCount, mTotalResiduals/mTotalCount, lowLLProbSampleCount, lowLLProbSampleCount/(double)mSampleCount));

    }

    private void logResidualData() {

        // report any stand-out sample or bucket, defined as X times the average
        double residualFactor = 3;

        double totalResidualsPerc = mTotalResiduals / mTotalCount;
        double resPercThreshold = totalResidualsPerc * residualFactor;
        int sampleCountAboveThreshold = 0;
        double countsAboveThreshold = 0;

        for (int i = 0; i < mSampleCount; ++i) {

            double resCount = mSampleResiduals[i];
            double sampleCount = mSampleTotals[i];

            if (resCount == 0)
                continue;

            double resPerc = resCount / sampleCount;

            if (resPerc >= resPercThreshold) {

//                LOGGER.debug(String.format("sample(%d) residuals(count=%.0f perc=%.2f) above avg(%.3f)",
//                        i, resCount, resPerc, totalResidualsPerc));

                ++sampleCountAboveThreshold;
                countsAboveThreshold += resCount;
            }
        }

        LOGGER.debug(String.format("samples above avg residuals: count(%d perc=%.3f) residuals(%.0f perc=%.3f))",
                sampleCountAboveThreshold, sampleCountAboveThreshold/(double)mSampleCount,
                countsAboveThreshold, countsAboveThreshold/mTotalResiduals));

        double bucketCountThreshold = mTotalCount / (double) mBucketCount * 0.25;

        for (int i = 0; i < mBucketCount; ++i) {

            double resCount = mBucketResiduals[i];
            double bucketCount = mBucketTotals[i];

            if (resCount == 0)
                continue;

            double resPerc = resCount / bucketCount;

            if (resPerc >= resPercThreshold && bucketCount >= bucketCountThreshold) {

                LOGGER.debug(String.format("bucket(%d) residuals(count=%.0f perc=%.2f) above avg(%.3f)",
                        i, resCount, resPerc, totalResidualsPerc));
            }
        }
    }

    private static boolean poissonProbEarlyExit(int mean, int observed)
    {
        if(mean <= 10)
            return (observed < 20);

        if(mean <= 100)
            return (observed >= mean - 30 && observed <= mean + 30);

        if(mean <= 1000)
            return (observed >= mean - 100 && observed <= mean + 100);

        if(mean <= 10000)
            return (observed >= mean - 300 && observed <= mean + 300);

        if(mean <= 100000)
            return (observed >= mean - 1000 && observed <= mean + 1000);

        if(mean <= 1000000)
            return (observed >= mean - 3000 && observed <= mean + 3000);

        return (observed >= mean - 5000 && observed <= mean + 5000);
    }


}
