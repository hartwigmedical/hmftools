package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.log;

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
    double[] mBucketTotals; // to help with seeding

    private NmfMatrix mW; // the bucket-signature values (x=BucketCount, y=SigCount)
    private NmfMatrix mH; // the sample-signature contributions (x=SigCount, y=SampleCount)
    private NmfMatrix mV; // the fitted matrix of samples and bucket counts (W x H)

    private NmfMatrix mRefSignatures;
    private NmfMatrix mRefContributions;
    double mSigFloatRate;
    private List<NmfMatrix> mStartSigs;

    // calculated values
    private double mTotalResiduals;
    private double mTotalCounts;
    private double[] mSampleTotals;
    private double[] mSampleResiduals;
    private boolean mIsValid;

    private Random mRandom;

    // internal constants
    private static double MIN_COST_CHANGE_PERCENT = 0.0001;

    public NmfCalculator(final NmfMatrix sampleBucketCounts, final NmfConfig config)
    {
        mConfig = config;

        mSigCount = mConfig.SigCount;
        mSampleCounts = sampleBucketCounts;
        mTotalCounts = mSampleCounts.sum();

        mBucketCount = sampleBucketCounts.Rows;
        mSampleCount = sampleBucketCounts.Cols;

        LOGGER.info("initialised with samples({}) buckets({}) sigCount({}) totalCount({})",
                mSampleCount, mBucketCount, mSigCount, mTotalCounts);

        mBucketTotals = new double[mBucketCount];

        for(int i = 0; i < mSampleCounts.Rows; ++i)
        {
            mBucketTotals[i] = DataUtils.sumVector(mSampleCounts.getRow(i));
        }

        mSampleTotals = new double[mSampleCounts.Cols];

        for(int i = 0; i < mSampleCounts.Cols; ++i)
        {
            mSampleTotals[i] = DataUtils.sumVector(mSampleCounts.getCol(i));
        }

        mTotalResiduals = 0;

        mW = null;
        mH = null;
        mV = new NmfMatrix(mBucketCount, mSampleCount);

        mRefSignatures = null;
        mRefContributions = null;
        mSigFloatRate = 1;
        mStartSigs = Lists.newArrayList();

        mSampleResiduals = new double[mSampleCount];
        mIsValid = false;

        // could seed from config instead
        mRandom = new Random(123456);
    }

    public void setSignatures(final NmfMatrix refSigs, double sigFloatRate)
    {
        mRefSignatures = refSigs;
        mSigFloatRate = sigFloatRate;
    }

    public void setContributions(final NmfMatrix refContributions) { mRefContributions = refContributions; }

    public final NmfMatrix getSignatures() { return mW; }
    public final NmfMatrix getContributions() { return mH; }
    public final NmfMatrix getFit() { return mV; }
    public double getTotalResiduals() { return mTotalResiduals; }
    public double getTotalCounts() { return mTotalCounts; }
    public List<NmfMatrix> getStartSigs() { return mStartSigs; }
    public final NmfMatrix getRefSignatures() { return mRefSignatures; }

    public boolean isValid() { return mIsValid; }

    public void performRun()
    {
        mIsValid = false;
        mTotalResiduals = 0;

        initSignatures();
        initContributions();
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

        // for each signature, either set a random bucket ratio or take the ref once
        // and add up all ratios, so they can be used to split the actual bucket count amongst them
        for (int s = 0; s < mSigCount; ++s) {

            bucketRatioTotal = 0;

            for (int i = 0; i < mW.Rows; ++i) {

                if (mRefContributions != null && s < mRefSignatures.Cols)
                    bucketRatios[i] = mRefSignatures.get(i, s);
                else
                    bucketRatios[i] = mRandom.nextDouble();

                bucketRatioTotal += bucketRatios[i];
            }

            for (int i = 0; i < mW.Rows; ++i) {
                wData[i][s] = bucketRatios[i] / bucketRatioTotal;
            }
        }

        mStartSigs.add(mW);
    }

    private void initContributions()
    {
        if(mRefContributions != null && mRefContributions.Rows == mSigCount)
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
        double sigTotal = 0;

        for (int n = 0; n < mSampleCounts.Cols; ++n) {

            double sampleTotal = mSampleTotals[n];
            sigTotal = 0;

            for (int s = 0; s < mSigCount; ++s) {

                if (mRefContributions != null && s < mRefContributions.Rows)
                    sigFractions[s] = mRefContributions.get(s, n);
                else
                    sigFractions[s] = mRandom.nextDouble();

                sigTotal += sigFractions[s];
            }

            for (int s = 0; s < mSigCount; ++s) {
                hData[s][n] = sampleTotal * sigFractions[s] / sigTotal;
            }
        }
    }

    private void calculate()
    {
        double prevCost = 0;
        double prevResiduals = 0;
        double prevDivergence = 0;

        int i = 0;
        for(; i < mConfig.MaxIterations; i++) {

            // compute output
            mV = mW.multiply(mH);

            // compare the original counts to the calculated matrix
            double cost = mSampleCounts.sumDiffSq(mV);

            if(Double.isNaN(cost))
            {
                LOGGER.warn("invalid cost value, exiting");
                mIsValid = false;
                break;
            }

            double costChange = prevCost > 0 ? (cost - prevCost) / prevCost : 0;

            if (mConfig.LogVerbose) {

                calcResiduals();
                double residualsChange = prevResiduals > 0 ? (mTotalResiduals - prevResiduals) / prevResiduals : 0;

                double divergence = calcDivergenceCost(false);
                double divergenceChange = prevDivergence > 0 ? (divergence - prevDivergence) / prevDivergence : 0;

                LOGGER.debug(String.format("%d: cost(%.0f chg=%.3f) residuals(%.1f chg=%.3f) divergence(%.0f chg=%.3f)",
                        i, cost, costChange, mTotalResiduals, residualsChange, divergence, divergenceChange));

                prevDivergence = divergence;
                prevResiduals = mTotalResiduals;
            }

            // check conditions to exit the fit routine
            if(i > 0)
            {
                if (cost < mConfig.ExitLevel)
                {
                    LOGGER.debug(String.format("%d: cost(%.1f) below cutoff(%.1f), exiting fit", i, cost, mConfig.ExitLevel));
                    break;

                } else if (abs(costChange) < MIN_COST_CHANGE_PERCENT) {

                    LOGGER.debug(String.format("%d: cost(%.1f -> %.1f), negligible change, exiting fit", i, cost, prevCost));
                    break;
                }
            }

            prevCost = cost;
            applyAdjustments();
        }

        if(i == mConfig.MaxIterations)
        {
            LOGGER.debug(String.format("%d: max iterations reached with cost(%.1f), exiting fit", i, prevCost));
        }

        calcResiduals();

        LOGGER.info(String.format("totalResiduals(%.1f) vs total(%.0f) as percent(%.3f)",
                mTotalResiduals, mTotalCounts, mTotalResiduals / mTotalCounts));

        mIsValid = true;
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
            double sigTotal = DataUtils.sumVector(mW.getCol(j));

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
            double sigTotal = DataUtils.sumVector(mH.getRow(j));

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

        final double[][] vData = mV.getData();

        for(int n = 0; n < mSampleCount; ++n)
        {
            double sampleResiduals = 0;

            for(int b = 0; b < mBucketCount; ++b)
            {
                double bucketCount = mSampleCounts.get(b, n);

                double sbContrib = vData[b][n];
                double absDiff = abs(bucketCount - sbContrib);
                sampleResiduals += absDiff;
            }

            mSampleResiduals[n] = sampleResiduals;
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

    private void logSampleData()
    {
        for(int n = 0; n < mSampleCount; ++n) {

            double sampleTotal = DataUtils.sumVector(mSampleCounts.getCol(n));

            LOGGER.debug(String.format("sample(%d) residuals(%.1f vs %.0f) total(%.1f)",
                    n, mSampleResiduals[n], sampleTotal, mTotalResiduals));
        }

    }
}
