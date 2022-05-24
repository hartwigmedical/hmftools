package com.hartwig.hmftools.sigs.nmf;

import static java.lang.Double.max;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.MatrixUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfCalculator
{

    private NmfConfig mConfig;

    private int mBucketCount;
    private int mSampleCount;
    private int mSigCount;

    private static final Logger LOGGER = LogManager.getLogger(NmfCalculator.class);

    // primary input - bucket counts per sample
    final Matrix mSampleCounts;
    private double mTotalCount;
    private double[] mBucketTotals; // to help with seeding
    private double[] mSampleTotals;

    private int mRunId;
    private Matrix mW; // the bucket-signature values (x=BucketCount, y=SigCount)
    private Matrix mH; // the sample-signature contributions (x=SigCount, y=SampleCount)
    private Matrix mV; // the fitted matrix of samples and bucket counts (W x H)
    private Matrix mPrevW;
    private Matrix mPrevH;
    private Matrix mPrevV;
    private boolean mIsValid;

    private Matrix mRefSignatures;
    private Matrix mRefContributions;
    private List<Matrix> mStartSigs;
    private Matrix mRandomStartSignatures;

    // calculated values
    private double mTotalResiduals;
    private double mNetResiduals;
    private double mLowestCost;

    private Random mRandom;

    // internal constants
    private static double MIN_COST_CHANGE_PERCENT = 0.00001;

    public NmfCalculator(final Matrix sampleBucketCounts, final NmfConfig config)
    {
        mConfig = config;
        mRunId = 0;

        mSigCount = 0; // will be set for each run
        mSampleCounts = sampleBucketCounts;
        mTotalCount = MatrixUtils.sum(mSampleCounts);

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
        mV = new Matrix(mBucketCount, mSampleCount);
        mPrevV = new Matrix(mBucketCount, mSampleCount);

        mRefSignatures = null;
        mRefContributions = null;
        mRandomStartSignatures = null;
        mStartSigs = Lists.newArrayList();

        mIsValid = false;

        // could seed from config instead
        mRandom = new Random(123456);
    }

    public void setSigCount(int sigCount) { mSigCount = sigCount; }

    public void setSignatures(final Matrix refSigs)
    {
        mRefSignatures = refSigs;
    }

    public void setContributions(final Matrix refContributions) { mRefContributions = refContributions; }

    public void setRandomSignatures(final Matrix randomSigs) { mRandomStartSignatures = randomSigs; }

    public final Matrix getSignatures() { return mW; }
    public final Matrix getContributions() { return mH; }
    public final Matrix getFit() { return mV; }
    public final Matrix getSampleCounts() { return mSampleCounts; }
    public double[] getBucketTotals() { return mBucketTotals; }
    public double[] getSampleTotals() { return mSampleTotals; }
    public double getTotalResiduals() { return mTotalResiduals; }
    public void clearLowestCost() { mLowestCost = 0; }

    public double getTotalCount() { return mTotalCount; }
    public final Matrix getRefSignatures() { return mRefSignatures; }

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
            MatrixUtils.multiply(mW, mH, mV, true);
            calcResiduals();

            LOGGER.debug(String.format("run %d: pre-fit: totalResiduals(%.0f) vs total(%.0f) as percent(%.5f)",
                    mRunId, mTotalResiduals, mTotalCount, mTotalResiduals / mTotalCount));
        }

        mPrevW = new Matrix(mBucketCount, mSigCount);
        mPrevH = new Matrix(mSigCount, mSampleCount);

        calculate();
    }

    private void initSignatures()
    {
        // to stick with convention, the signatures matrix is comprised of values between 0 - 1, with a sig's bucket ratios adding to 1
        // whereas the contributions per samples are its bucket counts split across the sigs
        if (mRefSignatures != null && mRefSignatures.Cols == mSigCount)
        {
            mW = new Matrix(mRefSignatures);
            return;
        }

        mW = new Matrix(mBucketCount, mSigCount);

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
            mH = new Matrix(mRefContributions);
            return;
        }

        mH = new Matrix(mSigCount, mSampleCount);

        // if the signatures are fractions of 1 for each bucket
        // then the contributions should be based around the actual bucket counts per sample
        // but divided randomly amongst the signatures
        double[][] hData = mH.getData();

        double[] sigFractions = new double[mSigCount];

        // if there are proposed or ref contributions in use, the other sigs should
        // be given a relatively small value compared to the ref
        // for now, assume there is only 1 proposed sig in play per sample
        double refSigAllocation = 1 - mConfig.SigFloatRate;
        // double refSigAllocation = 0.99;

        // non-proposed sigs need a contribution above zero to allow them to float
        double nonRefSigPercent = (1-refSigAllocation) / (mSigCount - 1);

        for (int n = 0; n < mSampleCounts.Cols; ++n)
        {
            double sampleTotal = mSampleTotals[n];
            double sigTotal = 0;

            // if this sample has specific contributions set, use the min contribution value for the restr
            // if not, set them all to randoms
            boolean requiresRandoms = true;

            if (mRefContributions != null)
            {
                for (int s = 0; s < mRefContributions.Rows; ++s)
                {
                    if (mRefContributions.get(s, n) > 0)
                    {
                        requiresRandoms = false;
                        break;
                    }
                }
            }

            for (int s = 0; s < mSigCount; ++s)
            {
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

            for (int s = 0; s < mSigCount; ++s)
            {
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

        for(; i < maxIterations; i++)
        {
            // compute the fit
            produceFit();

            if(mConfig.LogVerbose && i > 0)
            {
                logMatrixDiffs();
            }

            // compare the original counts to the calculated matrix
            currentCost = MatrixUtils.sumDiffSq(mSampleCounts, mV);

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

            if (mConfig.LogVerbose)
            {
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

                }
                else if (abs(costChange) < MIN_COST_CHANGE_PERCENT || costChange > 0)
                {

                    LOGGER.debug(String.format("run=%d, it=%d: cost(%.0f -> %.0f) percent(%.4f), %s change, exiting fit",
                            mRunId, i, prevCost, currentCost, costChange, costChange > 0 ? "positive" : "small"));
                    break;
                }

                // also check the rate of change to project whether it is likely to reach the current lowest cost level
                if(i > 10 && mLowestCost > 0)
                {
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

                    if (targetCostReduced > projectCostExit || targetCostLinear > projectCostExit)
                    {
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

        produceFit(); // ensure fit is the latest
        normaliseSignatures();

        if(!mIsValid)
            return;

        calcResiduals();

        mLowestCost = mLowestCost == 0 ? currentCost : min(mLowestCost, currentCost);

        LOGGER.info(String.format("run=%d, it=%d: residuals(%.0f) vs total(%.0f) as percent(%.5f) cost(init=%.0f early=%.0f end=%.0f lastChg=%.5f)",
                mRunId, i, mTotalResiduals, mTotalCount, mTotalResiduals / mTotalCount,
                initCost, earlyCost, currentCost, prevCostChange));
    }

    public void produceFit()
    {
        MatrixUtils.multiply(mW, mH, mV, true); // ensure fit is the latest
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
        Matrix wt = mW.transpose();
        Matrix hAdj = MatrixUtils.multiply(wt, mSampleCounts);
        Matrix hd = MatrixUtils.multiply(wt, mV);

        MatrixUtils.scalarDivide(hAdj, hd, true);
        MatrixUtils.scalarMultiply(mH, hAdj);

        if(mConfig.SigFloatRate > 0)
        {
            // update signatures matrix
            Matrix ht = mH.transpose();
            Matrix wAdj = MatrixUtils.multiply(mSampleCounts, ht);
            Matrix wd1 = MatrixUtils.multiply(mW, mH);
            Matrix wd = MatrixUtils.multiply(wd1, ht);

            MatrixUtils.scalarDivide(wAdj, wd, true);

            if(mConfig.SigFloatRate == 1)
            {
                MatrixUtils.scalarMultiply(mW, wAdj);
            }
            else
            {
                MatrixUtils.scalarMultiplyRateAdjusted(mW, wAdj, mConfig.SigFloatRate, mRefSignatures.Cols);
            }
        }
    }

    private void modelBrunet()
    {
        Matrix vWH = mSampleCounts;
        MatrixUtils.scalarDivide(vWH, mV);

        Matrix wSum = new Matrix(mSigCount, mSampleCount);
        double[][] wSumData = wSum.getData();
        for(int j = 0; j < mSigCount; ++j)
        {
            double sigTotal = sumVector(mW.getCol(j));

            for(int k = 0; k < mSampleCount; ++k)
            {
                wSumData[j][k] = sigTotal;
            }
        }

        Matrix wt = mW.transpose();
        Matrix wt_vWH = MatrixUtils.multiply(wt, vWH);
        Matrix hAdj = wt_vWH;
        MatrixUtils.scalarDivide(hAdj, wSum);

        MatrixUtils.scalarMultiply(mH, hAdj);

        // recalc V and WH using the new H
        mV = MatrixUtils.multiply(mW, mH);
        vWH = mSampleCounts;
        MatrixUtils.scalarDivide(vWH, mV);

        // now adjust W
        Matrix hSum = new Matrix(mBucketCount, mSigCount);
        double[][] hSumData = hSum.getData();
        for(int j = 0; j < mSigCount; ++j)
        {
            double sigTotal = sumVector(mH.getRow(j));

            for(int k = 0; k < mBucketCount; ++k)
            {
                hSumData[k][j] = sigTotal;
            }
        }

        Matrix ht = mH.transpose();
        Matrix vWH_ht = MatrixUtils.multiply(vWH, ht);
        Matrix wAdj = vWH_ht;
        MatrixUtils.scalarDivide(wAdj, hSum);

        MatrixUtils.scalarMultiply(mW, wAdj);
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

        for (int n = 0; n < mSampleCount; ++n)
        {
            for (int b = 0; b < mBucketCount; ++b)
            {
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

            for(int j = 0; j < mH.Cols; ++j)
            {
                hData[i][j] *= contribAdj;
            }

            // verify bucket ratios for a sig now total 1 and so are in effect percentages
            if(!doublesEqual(sumVector(mW.getCol(i)),1))
            {
                LOGGER.error("sig normalisation failed");
                mIsValid = false;
                return;
            }
        }

        // check that V hasn't changed
        final Matrix vCopy = new Matrix(mV);
        produceFit();

        if(!mV.hasValidData(false))
        {
            LOGGER.warn("V-matrix has invalid data");
            mIsValid = false;
            return;
        }

        double sumDiff = MatrixUtils.sumDiffSq(mV, vCopy);
        boolean matrixEqual = mV.equals(vCopy);

        if(!doublesEqual(sumDiff, 0) || !matrixEqual)
        {
            LOGGER.warn(String.format("bucket ratio adjustments failed: diff(%.4f) equal(%s)",
                    sumDiff, matrixEqual));
            mIsValid = false;
        }
    }

    private void logMatrixDiffs()
    {
        Matrix relDiff = MatrixUtils.getDiff(mV, mPrevV, true);
        Matrix absDiff = MatrixUtils.getDiff(mV, mPrevV, false);
        double avgPercChange = MatrixUtils.sum(relDiff) / (mV.Rows * mV.Cols);
        LOGGER.debug(String.format("V-matrix diffs: abs(%.0f) relative(%.4f)", MatrixUtils.sum(absDiff), avgPercChange));

        relDiff = MatrixUtils.getDiff(mW, mPrevW, true);
        absDiff = MatrixUtils.getDiff(mW, mPrevW, false);
        avgPercChange = MatrixUtils.sum(relDiff) / (mW.Rows * mW.Cols);
        LOGGER.debug(String.format("W-matrix diffs: abs(%.0f) relative(%.4f)", MatrixUtils.sum(absDiff), avgPercChange));

        relDiff = MatrixUtils.getDiff(mH, mPrevH, true);
        absDiff = MatrixUtils.getDiff(mH, mPrevH, false);
        avgPercChange = MatrixUtils.sum(relDiff) / (mH.Rows * mH.Cols);
        LOGGER.debug(String.format("H-matrix diffs: abs(%.0f) relative(%.4f)", MatrixUtils.sum(absDiff), avgPercChange));

    }
}
