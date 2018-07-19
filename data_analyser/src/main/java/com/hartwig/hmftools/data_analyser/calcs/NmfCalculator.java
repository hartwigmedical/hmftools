package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Double.max;
import static java.lang.Double.valueOf;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.min;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfCalculator {

    /*
     * matrix input - bucket counts (b) x samples (n)
     * w.mult(h) = this
     * W is weights matrix - bucket (b) x sigs (r)
     * H is feature matrix - sigs (r) x samples (n)
     *
     * r = sig count (size of factor matrices)
     * maxIterations is NUMBER OF ITERATIONS.
     * e is ERROR (once we're below that, we can return)
     *
     * returns W and H matrices
     */

    private NmfConfig mConfig;

    private int mBucketCount;
    private int mSampleCount;
    private int mSigCount;

    private static final Logger LOGGER = LogManager.getLogger(NmfCalculator.class);

    // primary input - bucket counts per sample
    final NmfMatrix mSampleCounts;
    private double mTotalCount;
    private double[] mBucketTotals; // to help with seeding
    private double[] mSampleTotals;

    private int mRunId;
    private NmfMatrix mW; // the bucket-signature values (x=BucketCount, y=SigCount)
    private NmfMatrix mH; // the sample-signature contributions (x=SigCount, y=SampleCount)
    private NmfMatrix mV; // the fitted matrix of samples and bucket counts (W x H)
    private NmfMatrix mPrevW;
    private NmfMatrix mPrevH;
    private NmfMatrix mPrevV;
    private boolean mIsValid;

    private NmfMatrix mRefSignatures;
    private NmfMatrix mRefContributions;
    private double mSigFloatRate;
    private List<NmfMatrix> mStartSigs;
    private NmfMatrix mRandomStartSignatures;

    // calculated values
    private double mTotalResiduals;
    private double mNetResiduals;
    private double mLowestCost;

    private Random mRandom;

    // internal constants
    private static double MIN_COST_CHANGE_PERCENT = 0.00001;

    public NmfCalculator(final NmfMatrix sampleBucketCounts, final NmfConfig config)
    {
        mConfig = config;
        mRunId = 0;

        mSigCount = 0; // will be set for each run
        mSampleCounts = sampleBucketCounts;
        mTotalCount = mSampleCounts.sum();

        mBucketCount = sampleBucketCounts.Rows;
        mSampleCount = sampleBucketCounts.Cols;

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

        mTotalResiduals = 0;
        mNetResiduals = 0;
        mLowestCost = 0;

        mW = null;
        mH = null;
        mV = new NmfMatrix(mBucketCount, mSampleCount);
        mPrevV = new NmfMatrix(mBucketCount, mSampleCount);

        mRefSignatures = null;
        mRefContributions = null;
        mRandomStartSignatures = null;
        mSigFloatRate = 1;
        mStartSigs = Lists.newArrayList();

        mIsValid = false;

        // could seed from config instead
        mRandom = new Random(123456);
    }

    public void setSigCount(int sigCount) { mSigCount = sigCount; }

    public void setSignatures(final NmfMatrix refSigs, double sigFloatRate)
    {
        mRefSignatures = refSigs;
        mSigFloatRate = sigFloatRate;
    }

    public void setContributions(final NmfMatrix refContributions) { mRefContributions = refContributions; }

    public void setRandomSignatures(final NmfMatrix randomSigs) { mRandomStartSignatures = randomSigs; }

    public final NmfMatrix getSignatures() { return mW; }
    public final NmfMatrix getContributions() { return mH; }
    public final NmfMatrix getFit() { return mV; }
    public final NmfMatrix getSampleCounts() { return mSampleCounts; }
    public double[] getBucketTotals() { return mBucketTotals; }
    public double[] getSampleTotals() { return mSampleTotals; }
    public double getTotalResiduals() { return mTotalResiduals; }
    public void clearLowestCost() { mLowestCost = 0; }

    public double getTotalCount() { return mTotalCount; }
    public List<NmfMatrix> getStartSigs() { return mStartSigs; }
    public final NmfMatrix getRefSignatures() { return mRefSignatures; }

    public boolean isValid() { return mIsValid; }

    public void performRun(final int runId)
    {
        mRunId = runId;
        mIsValid = false;
        mTotalResiduals = 0;

        if(mSigCount <= 0)
            return;

//        LOGGER.debug("initialised with samples({}) buckets({}) sigCount({}) totalCount({})",
//                mSampleCount, mBucketCount, mSigCount, mTotalCount);

        initSignatures();
        initContributions();

        if(mConfig.LogVerbose && mRefSignatures != null && mRefContributions != null)
        {
            mW.multiply(mH, mV, true);
            calcResiduals();

            LOGGER.debug(String.format("run %d: pre-fit: totalResiduals(%.0f) vs total(%.0f) as percent(%.5f)",
                    mRunId, mTotalResiduals, mTotalCount, mTotalResiduals / mTotalCount));
        }

        mPrevW = new NmfMatrix(mBucketCount, mSigCount);
        mPrevH = new NmfMatrix(mSigCount, mSampleCount);

        calculate();
    }

    private void initSignatures()
    {
        // to stick with convention, the signatures matrix is comprised of values between 0 - 1, with a sig's bucket ratios adding to 1
        // whereas the contributions per samples are its bucket counts split across the sigs
        if (mRefSignatures != null && mRefSignatures.Cols == mSigCount) {

            mW = new NmfMatrix(mRefSignatures);
            return;
        }

        mW = new NmfMatrix(mBucketCount, mSigCount);

        double[][] wData = mW.getData();

        double[] bucketRatios = new double[mBucketCount];
        double bucketRatioTotal = 0;

        List<Integer> randomSigs = Lists.newArrayList();
        if(mRandomStartSignatures != null)
        {
            // get the next set of random sigs
            int nextSigIndex = mRandom.nextInt(mRandomStartSignatures.Cols);
            int attempts = 0;
            while(randomSigs.size() < mSigCount && attempts < mSigCount * 5)
            {
                if(!randomSigs.contains(nextSigIndex))
                {
                    randomSigs.add(nextSigIndex);
                }

                nextSigIndex = mRandom.nextInt(mRandomStartSignatures.Cols);
                ++attempts;
            }

            if(randomSigs.size() < mSigCount)
            {
                LOGGER.warn("insufficient random sigs found");
                return;
            }

            // LOGGER.debug("using random sigs: {}", randomSigs.toString());
        }

        // for each signature, either set a random bucket ratio or take the ref one
        // and add up all ratios, so they can be used to split the actual bucket count amongst them
        int nonRefSigIndex = 0;
        for (int s = 0; s < mSigCount; ++s) {

            bucketRatioTotal = 0;

            if (mRefContributions != null && s < mRefSignatures.Cols) {

                for (int i = 0; i < mW.Rows; ++i) {
                    bucketRatios[i] = mRefSignatures.get(i, s);
                    bucketRatioTotal += bucketRatios[i];
                }
            }
            else if(mRandomStartSignatures != null)
            {
                int randomSig = randomSigs.get(nonRefSigIndex);

                for (int i = 0; i < mW.Rows; ++i) {

                    bucketRatios[i] = mRandomStartSignatures.get(i, randomSig);
                    bucketRatioTotal += bucketRatios[i];
                }

                ++nonRefSigIndex;
            }
            else
            {
                // random values
                for (int i = 0; i < mW.Rows; ++i) {

                    bucketRatios[i] = mRandom.nextDouble();
                    bucketRatioTotal += bucketRatios[i];
                }
            }

            // ensure sig bucket ratios total to 1 (by convention)
            for (int i = 0; i < mW.Rows; ++i) {
                wData[i][s] = bucketRatios[i] / bucketRatioTotal;
            }
        }

        // mStartSigs.add(mW);
    }

    private void initContributions()
    {
        if(mRefContributions != null && mRefContributions.Rows == mSigCount && mConfig.UseRefSigs)
        {
            mH = new NmfMatrix(mRefContributions);
            return;
        }

        mH = new NmfMatrix(mSigCount, mSampleCount);

        // if the signatures are fractions of 1 for each bucket
        // then the contributions should be based around the actual bucket counts per sample
        // but divided randomly amongst the signatures
        double[][] hData = mH.getData();

        double[] sigFractions = new double[mSigCount];

        // if there are proposed or ref contributions in use, the other sigs should
        // be given a relatively small value compared to the ref
        // for now, assume there is only 1 proposed sig in play per sample
        double refSigAllocation = 1 - mSigFloatRate;
        // double refSigAllocation = 0.99;

        // non-proposed sigs need a contribution above zero to allow them to float
        double nonRefSigPercent = (1-refSigAllocation) / (mSigCount - 1);

        for (int n = 0; n < mSampleCounts.Cols; ++n) {

            double sampleTotal = mSampleTotals[n];
            double sigTotal = 0;

            // if this sample has specific contributions set, use the min contribution value for the restr
            // if not, set them all to randoms
            boolean requiresRandoms = true;

            if (mRefContributions != null) {

                for (int s = 0; s < mRefContributions.Rows; ++s) {

                    if (mRefContributions.get(s, n) > 0) {

                        requiresRandoms = false;
                        break;
                    }
                }
            }

            for (int s = 0; s < mSigCount; ++s) {

                if(requiresRandoms)
                {
                    sigFractions[s] = mRandom.nextDouble();
                }
                else
                {
                    // each sample is given either designated sig contribution or the suitable low value calc'ed above
                    if (s < mRefContributions.Rows)
                    {
                        sigFractions[s] = max(mRefContributions.get(s, n), nonRefSigPercent);
                    }
                    else
                    {
                        sigFractions[s] = nonRefSigPercent;
                    }
                }

                sigTotal += sigFractions[s];
            }

            for (int s = 0; s < mSigCount; ++s) {
                hData[s][n] = sampleTotal * sigFractions[s] / sigTotal;
            }
        }
    }

    private void calculate()
    {
        mIsValid = true;

        double currentCost = 0;
        double prevCost = 0;
        double prevResiduals = 0;
        double prevDivergence = 0;
        double prevCostChange = 0;
        double initCost = 0;
        double earlyCost = 0;

        int i = 0;
        int iterCheckInterval = 10; // how often to check, rather than every time
        int maxIterations = mConfig.MaxIterations;
        int permittedExtensions = 3;
        double projectCostExit = mLowestCost * 1.25; // build a buffer in for uncertainty

        for(; i < maxIterations; i++) {

            // compute the fit
            mW.multiply(mH, mV, true);

            if(mConfig.LogVerbose && i > 0) {
                logMatrixDiffs();
            }

            // compare the original counts to the calculated matrix
            currentCost = mSampleCounts.sumDiffSq(mV);

            if(i == 0)
                initCost = currentCost;
            else if(i == iterCheckInterval)
                earlyCost = currentCost;

            if(Double.isNaN(currentCost) || Double.isInfinite(currentCost) || currentCost > 1e50)
            {
                LOGGER.warn("it={}: invalid cost value: nan={} infinite={} max={}, exiting",
                        i, Double.isNaN(currentCost), Double.isInfinite(currentCost), currentCost > 1e50);

                mIsValid = false;
                break;
            }

            double costChange = prevCost > 0 ? (currentCost - prevCost) / prevCost : 0;

            if (mConfig.LogVerbose) {

                calcResiduals();
                double residualsChange = prevResiduals > 0 ? (mTotalResiduals - prevResiduals) / prevResiduals : 0;

                double divergence = calcDivergenceCost(false);
                double divergenceChange = prevDivergence > 0 ? (divergence - prevDivergence) / prevDivergence : 0;

                LOGGER.debug(String.format("%d: cost(%.0f chg=%.3f) residuals(%.1f chg=%.3f) divergence(%.0f chg=%.3f)",
                        i, currentCost, costChange, mTotalResiduals, residualsChange, divergence, divergenceChange));

                prevDivergence = divergence;
                prevResiduals = mTotalResiduals;
            }

            // check conditions to exit the fit routine
            if(i > 0 && (i % iterCheckInterval) == 0)
            {
                if (currentCost < mConfig.ExitLevel)
                {
                    LOGGER.debug(String.format("run=%d, it=%d: cost(%.0f) below cutoff(%.0f), exiting fit", mRunId, i, currentCost, mConfig.ExitLevel));
                    break;

                } else if (abs(costChange) < MIN_COST_CHANGE_PERCENT || costChange > 0) {

                    LOGGER.debug(String.format("run=%d, it=%d: cost(%.0f -> %.0f) percent(%.4f), %s change, exiting fit",
                            mRunId, i, prevCost, currentCost, costChange, costChange > 0 ? "positive" : "small"));
                    break;
                }

                // also check the rate of change to project whether it is likely to reach the current lowest cost level
                if(i > 10 && mLowestCost > 0) {

                    double changeRate = (prevCostChange - costChange) / prevCostChange;
                    int remainingIts = mConfig.MaxIterations - i;

                    // firstly optimistically assume the current reduction rate continues
                    double targetCostLinear = currentCost * Math.pow(1 + costChange, remainingIts);

                    double targetCostReduced = targetCostLinear;

                    if(changeRate > 0 && changeRate < 0.05)
                    {
                        // otherwise assume it also slows, so work out how many iterations until it reaches the flat line threshold
                        // but here even keeping the rate of change constant is optimistic
                        double minChangeIts = min(log(MIN_COST_CHANGE_PERCENT) / log(1 - changeRate), remainingIts);
                        targetCostReduced = currentCost * Math.pow(1 + costChange, minChangeIts);
                    }

                    if (targetCostReduced > projectCostExit || targetCostLinear > projectCostExit) {

                        LOGGER.debug(String.format(
                                "run=%d, it=%d: costChange(%.6f percVsLast=%.4f) to small for cost(%.0f vs low=%.0f) projected(lin=%.0f red=%.0f), exiting fit",
                                mRunId, i, costChange, changeRate, currentCost, mLowestCost, targetCostLinear, targetCostReduced));
                        break;
                    }
                }
            }

            prevCost = currentCost;
            prevCostChange = costChange;

            if(mConfig.LogVerbose) {

                mPrevV.setData(mV.getData());
                mPrevW.setData(mW.getData());
                mPrevH.setData(mH.getData());
            }

            applyAdjustments();

            if(i == maxIterations - 1)
            {
                // prior to exiting, check if worth continuing on if the current run is already the best fit
                if(mLowestCost > 0 && currentCost < mLowestCost && permittedExtensions > 0)
                {
                    LOGGER.debug(String.format("run=%d, it=%d: extending max iterations with new lowest cost(%.0f vs prev=%.0f) change(%.4f)",
                            mRunId, i, currentCost, mLowestCost, costChange));

                    maxIterations += mConfig.MaxIterations;
                    --permittedExtensions;
                }
                else
                {
                    LOGGER.debug(String.format("run=%d, it=%d: max iterations reached with cost(%.0f) change(%.6f), exiting fit",
                            mRunId, i, currentCost, costChange));
                    break;
                }
            }
        }

        if(!mIsValid || !mW.hasValidData(false) || !mH.hasValidData(false) || !mV.hasValidData(false))
            return;

        mW.multiply(mH, mV, true); // ensure fit is the latest
        normaliseSignatures();
        calcResiduals();

        mLowestCost = mLowestCost == 0 ? currentCost : min(mLowestCost, currentCost);

        LOGGER.info(String.format("run=%d, it=%d: residuals(%.0f) vs total(%.0f) as percent(%.5f) cost(init=%.0f early=%.0f end=%.0f lastChg=%.5f)",
                mRunId, i, mTotalResiduals, mTotalCount, mTotalResiduals / mTotalCount,
                initCost, earlyCost, currentCost, prevCostChange));
    }

    private void applyAdjustments()
    {
        switch(mConfig.Model)
        {
            case BRUNET:
                modelBrunet();
                break;

            case STANDARD:
            default:
                modelStandard();
                break;
        }
    }

    private void modelStandard()
    {
        // the multiplicative update method (described by Lee and Seund, 2001)
        // https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization.pdf

        // update contribution matrix
        NmfMatrix wt = mW.transpose();
        NmfMatrix hAdj = wt.multiply(mSampleCounts);
        NmfMatrix hd = wt.multiply(mV);

        hAdj.scalarDivide(hd);
        mH.scalarMultiply(hAdj);

        if(mSigFloatRate > 0) {

            // update signatures matrix
            NmfMatrix ht = mH.transpose();
            NmfMatrix wAdj = mSampleCounts.multiply(ht);
            NmfMatrix wd1 = mW.multiply(mH);
            NmfMatrix wd = wd1.multiply(ht);

            wAdj.scalarDivide(wd);

            if(mSigFloatRate == 1)
            {
                mW.scalarMultiply(wAdj);
            }
            else
            {
                mW.scalarMultiplyRateAdjusted(wAdj, mSigFloatRate, mRefSignatures.Cols);
            }
        }
    }

    private void modelBrunet()
    {
        NmfMatrix vWH = mSampleCounts;
        vWH.scalarDivide(mV);

        NmfMatrix wSum = new NmfMatrix(mSigCount, mSampleCount);
        double[][] wSumData = wSum.getData();
        for(int j = 0; j < mSigCount; ++j)
        {
            double sigTotal = sumVector(mW.getCol(j));

            for(int k = 0; k < mSampleCount; ++k)
            {
                wSumData[j][k] = sigTotal;
            }
        }

        NmfMatrix wt = mW.transpose();
        NmfMatrix wt_vWH = wt.multiply(vWH);
        NmfMatrix hAdj = wt_vWH;
        hAdj.scalarDivide(wSum);

        mH.scalarMultiply(hAdj);

        // recalc V and WH using the new H
        mV = mW.multiply(mH);
        vWH = mSampleCounts;
        vWH.scalarDivide(mV);

        // now adjust W
        NmfMatrix hSum = new NmfMatrix(mBucketCount, mSigCount);
        double[][] hSumData = hSum.getData();
        for(int j = 0; j < mSigCount; ++j)
        {
            double sigTotal = sumVector(mH.getRow(j));

            for(int k = 0; k < mBucketCount; ++k)
            {
                hSumData[k][j] = sigTotal;
            }
        }

        NmfMatrix ht = mH.transpose();
        NmfMatrix vWH_ht = vWH.multiply(ht);
        NmfMatrix wAdj = vWH_ht;
        wAdj.scalarDivide(hSum);

        mW.scalarMultiply(wAdj);
    }

    private void calcResiduals()
    {
        mTotalResiduals = 0;
        mNetResiduals = 0;

        final double[][] vData = mV.getData();
        final double[][] scData = mSampleCounts.getData();

        for(int n = 0; n < mSampleCount; ++n)
        {
            double sampleResiduals = 0;

            for(int b = 0; b < mBucketCount; ++b)
            {
                double bucketCount = scData[b][n];

                double sbContrib = vData[b][n];
                double diff = bucketCount - sbContrib;
                double absDiff = abs(diff);
                sampleResiduals += absDiff;
                mNetResiduals += diff;
            }

            mTotalResiduals += sampleResiduals;
        }
    }


    private double calcDivergenceCost(boolean useVAsRef)
    {
        // Kullback-Leibler divergence: Aij * log(Aij/Bij) - Aij + Bij
        double divergSum = 0;

        final double[][] vData = mV.getData();
        final double[][] sbData = mSampleCounts.getData();

        for (int n = 0; n < mSampleCount; ++n) {

            for (int b = 0; b < mBucketCount; ++b) {

                double A = useVAsRef ? vData[b][n] : sbData[b][n];
                double B = !useVAsRef ? vData[b][n] : sbData[b][n];

                if(B == 0)
                    B = 0.001;

                if(A == 0)
                    A = 0.001;

                double diverg = A * log(A/B) - A + B;
                divergSum += diverg;
            }
        }

        return divergSum;
    }

    private void normaliseSignatures()
    {
        if(mConfig.SigFloatRate == 0)
            return;

        // adjust all signature bucket ratios to sum to 1, and adjust contributions accordingly
        double[][] wData = mW.getData();
        double[][] hData = mH.getData();

        for(int i = 0; i < mW.Cols; ++i)
        {
            double bucketRatioTotal = sumVector(mW.getCol(i));

            if(bucketRatioTotal == 0)
                continue;

            double contribAdj = 0;

            for(int j = 0; j < mW.Rows; ++j)
            {
                // bucket ratio: x -> x/total to make a percentage
                if(j == 0) {
                    double prevVal = wData[j][i];
                    wData[j][i] /= bucketRatioTotal;
                    contribAdj = prevVal / wData[j][i];
                }
                else
                {
                    wData[j][i] /= bucketRatioTotal;
                }
            }

            for(int j = 0; j < mH.Cols; ++j) {

                hData[i][j] *= contribAdj;
            }

                // verify bucket ratios for a sig now total 1 and so are in effect percentages
            if(!doublesEqual(sumVector(mW.getCol(i)),1))
            {
                LOGGER.error("sig normalisation failed");
                return;
            }
        }

        // check that V hasn't changed
        final NmfMatrix vCopy = new NmfMatrix(mV);
        mW.multiply(mH, mV, true);

        double sumDiff = mV.sumDiffSq(vCopy);
        boolean matrixEqual = mV.equals(vCopy);

        if(!doublesEqual(sumDiff, 0) || !matrixEqual)
        {
            LOGGER.warn(String.format("bucket ratio adjustments failed: diff(%.4f) equal(%s)",
                    sumDiff, matrixEqual));
        }
    }

    private void logMatrixDiffs()
    {
        NmfMatrix relDiff = NmfMatrix.getDiff(mV, mPrevV, true);
        NmfMatrix absDiff = NmfMatrix.getDiff(mV, mPrevV, false);
        double avgPercChange = relDiff.sum() / (mV.Rows * mV.Cols);
        LOGGER.debug(String.format("V-matrix diffs: abs(%.0f) relative(%.4f)", absDiff.sum(), avgPercChange));

        relDiff = NmfMatrix.getDiff(mW, mPrevW, true);
        absDiff = NmfMatrix.getDiff(mW, mPrevW, false);
        avgPercChange = relDiff.sum() / (mW.Rows * mW.Cols);
        LOGGER.debug(String.format("W-matrix diffs: abs(%.0f) relative(%.4f)", absDiff.sum(), avgPercChange));

        relDiff = NmfMatrix.getDiff(mH, mPrevH, true);
        absDiff = NmfMatrix.getDiff(mH, mPrevH, false);
        avgPercChange = relDiff.sum() / (mH.Rows * mH.Cols);
        LOGGER.debug(String.format("H-matrix diffs: abs(%.0f) relative(%.4f)", absDiff.sum(), avgPercChange));

    }
}
