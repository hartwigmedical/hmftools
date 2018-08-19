package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.calcAbsDiffs;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.initVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVectors;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NmfSampleFitter {

    private static final Logger LOGGER = LogManager.getLogger(NmfSampleFitter.class);

    private NmfConfig mConfig;
    private NmfMatrix mRefSignatures;
    private NmfMatrix mSampleCounts;
    private NmfMatrix mAllContributions;

    private int[] mSigCountFrequency;

    private boolean mIsValid;

    private static double MAX_FIT_CSS_DIFF = 0.01;
    private static double MINOR_SIG_CONTRIBUTION_PERC = 0.001;
    private static double MIN_SIG_CONTRIBUTION_PERC = 0.01;


    public NmfSampleFitter(final NmfConfig config, final NmfMatrix sampleCounts, final NmfMatrix refSigs)
    {
        mConfig = config;
        mSampleCounts = sampleCounts;
        mAllContributions = new NmfMatrix(refSigs.Cols, sampleCounts.Cols);
        mSigCountFrequency = new int[refSigs.Cols];
        mRefSignatures = refSigs;
        mIsValid = true;
    }

    public final NmfMatrix getContributions() { return mAllContributions; }
    public boolean isValid() { return mIsValid; }

    public void fitSamples()
    {
        for(int i = 0; i < mSampleCounts.Cols; ++i)
        {
            boolean sampleValid = fitSample(i);

            if(!sampleValid)
            {
                mIsValid = false;
                break;
            }
        }

        // report frequency of how many sigs are used across the cohort
        for(int i = 0; i < mSigCountFrequency.length; ++i)
        {
            if(mSigCountFrequency[i] == 0)
                continue;

            LOGGER.debug("assigned sig count({}) with frequency({})", i, mSigCountFrequency[i]);
        }

        SigReporter sigReporter = new SigReporter(mSampleCounts, mRefSignatures, mAllContributions, mRefSignatures, mRefSignatures, mConfig);
        sigReporter.runAnalysis();
    }

    private boolean fitSample(final int sampleId)
    {
        int bucketCount = mSampleCounts.Rows;

        NmfMatrix sampleMatrix = new NmfMatrix(bucketCount, 1);

        final double[] sampleCounts = mSampleCounts.getCol(sampleId);
        double[][] sData = sampleMatrix.getData();
        double sampleCount = 0;

        for(int i = 0; i < bucketCount; ++i)
        {
            sData[i][0] = sampleCounts[i];
            sampleCount += sampleCounts[i];
        }

        int refSigCount = mRefSignatures.Cols;

        boolean[] sigsInUse = new boolean[refSigCount];
        NmfMatrix reducedSigs = new NmfMatrix(mRefSignatures);
        int currentSigCount = 0;

        // start with all in use unless below required ML threshold
        for(int i = 0; i < refSigCount; ++i)
        {
            if(belowRequiredMutLoad(i, sampleCount))
            {
                reduceSigData(sigsInUse, reducedSigs, i);
            }
            else
            {
                sigsInUse[i] = true;
                ++currentSigCount;
            }
        }

        NmfCalculator nmfCalc = new NmfCalculator(sampleMatrix, mConfig);

        nmfCalc.setSigCount(refSigCount);

        double[] prevContribs = new double[refSigCount];
        double prevResiduals = 0;
        double lastFitVsActualCss = 0;
        int lastSigRemoved = 0;
        int iterations = 0;

        while(currentSigCount >= 1)
        {
            nmfCalc.setSignatures(reducedSigs);
            nmfCalc.performRun(sampleId);

            if(!nmfCalc.isValid())
            {
                LOGGER.warn("NMF fit failed for sample({})", sampleId);
                return false;
            }

            final double[] newContribs = nmfCalc.getContributions().getCol(0);
            final double[] newFit = nmfCalc.getFit().getCol(0);

            double newFitVsActualCss = calcCSS(newFit, sampleCounts);

            if(iterations > 0) {

                // assess the new contributions vs the actual counts, and exit if the change is greater than 0.01
                double cssDiff = lastFitVsActualCss - newFitVsActualCss;

                if (cssDiff > MAX_FIT_CSS_DIFF) {

                    LOGGER.debug(String.format("sample(%d) fitVsActualsCss(%.4f -> %.4f diff=%.4f) with sigCount(%d), exiting",
                            sampleId, lastFitVsActualCss, newFitVsActualCss, cssDiff, currentSigCount));

                    // restore to last sig count prior to drop-off
                    ++currentSigCount;
                    enableRefSig(reducedSigs, lastSigRemoved);
                    nmfCalc.produceFit();
                    break;
                }
            }

            copyVector(newContribs, prevContribs);
            lastFitVsActualCss = newFitVsActualCss;
            prevResiduals = nmfCalc.getTotalResiduals();

            if(currentSigCount == 1)
            {
                LOGGER.debug(String.format("sample(%d) fitVsActualsCss(%.4f) single sig, exiting",
                        sampleId, lastFitVsActualCss));
                break;
            }

            // remove the least contributing signature
            int leastIndex = -1;
            double leastContrib = 0;
            int minorSigsRemoved = 0;

            for (int i = 0; i < refSigCount; ++i) {

                if(!sigsInUse[i])
                    continue;

                double sigContrib = newContribs[i];
                double sigPercent = sigContrib/sampleCount;

                if(sigContrib < 1 || sigPercent < MINOR_SIG_CONTRIBUTION_PERC)
                {
                    sigsInUse[i] = false;
                    ++minorSigsRemoved;

                    reduceSigData(sigsInUse, reducedSigs, i);
                    continue;
                }

                if(minorSigsRemoved == 0) { //  && sigPercent < MIN_SIG_CONTRIBUTION_PERC

                    // find the lowest contributing sig to remove (if no others have been this round already)
                    if (leastIndex == -1 || newContribs[i] < leastContrib) {
                        leastIndex = i;
                        leastContrib = prevContribs[i];
                    }
                }
            }

            if(minorSigsRemoved > 0)
            {
                currentSigCount -= minorSigsRemoved;

                LOGGER.debug(String.format("sample(%d) removed %d minor sigs, remaining(%d)", sampleId, minorSigsRemoved, currentSigCount));
            }
            else if(leastIndex >= 0)
            {
                --currentSigCount;
                lastSigRemoved = leastIndex;

                LOGGER.debug(String.format("sample(%d) remove least sig(%d contrib=%.1f perc=%.3f), fitVsActualsCss(%.4f) remaining(%d)",
                        sampleId, leastIndex, leastContrib, leastContrib/sampleCount, lastFitVsActualCss, currentSigCount));

                reduceSigData(sigsInUse, reducedSigs, leastIndex);
            }
            else
            {
                LOGGER.debug(String.format("sample(%d) no more sigs removed, count(%d)", sampleId, currentSigCount));
                break;
            }

            ++iterations;
        }

        int requiredSigsAdded = addRequiredSigs(sigsInUse, reducedSigs);

        if(requiredSigsAdded > 0)
        {
            currentSigCount += requiredSigsAdded;

            LOGGER.debug(String.format("sample(%d) final run with %d required sigs added back, sigCount(%d)",
                    sampleId, requiredSigsAdded, currentSigCount));

            nmfCalc.setSignatures(reducedSigs);
            nmfCalc.performRun(sampleId);
            prevResiduals = nmfCalc.getTotalResiduals();
            prevContribs = nmfCalc.getContributions().getCol(0);
        }

        // use the previous fit's contributions ie when still had a sufficiently high CSS
        double[][] allContribData = mAllContributions.getData();

        for(int i = 0; i < refSigCount; ++i)
        {
            allContribData[i][sampleId] = prevContribs[i];
        }

        double fitTotal = sumVector(nmfCalc.getFit().getCol(0));

        LOGGER.debug(String.format("sample(%d) residualsPerc(%.0f of %.0f fit=%.0f perc=%.4f) sigCount(%d)",
                sampleId, prevResiduals, sampleCount, fitTotal, prevResiduals/sampleCount, currentSigCount));

        if(currentSigCount < 1 || currentSigCount > refSigCount)
        {
            LOGGER.error("invalid sample sig count({})", currentSigCount);
        }
        else
        {
            mSigCountFrequency[currentSigCount] += 1;
        }

        return true;
    }

    private void reduceSigData(boolean[] sigsInUse, NmfMatrix sigs, int sigIndex)
    {
        sigsInUse[sigIndex] = false;

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex, 0);
        }
    }

    private boolean belowRequiredMutLoad(int sig, double sampleCount)
    {
        switch(sig)
        {
            case 12:
            case 13:
                return sampleCount < 1e4;

            case 5:
            case 17:
            case 18:
            case 24:
            case 25:
            case 30:
            case 48:
                return sampleCount < 1e4;

            default:
                return false;
        }

    }

    private void enableRefSig(NmfMatrix sigs, int sigIndex)
    {
        final double[][] refData = mRefSignatures.getData();

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex,  refData[i][sigIndex]);
        }
    }

    private int addRequiredSigs(boolean[] sigsInUse, NmfMatrix sigs)
    {
        // force some signatures to be included: currently 2 and 5
        int sigsAdded = 0;

        if(!sigsInUse[0])
        {
            sigsInUse[0] = true;
            enableRefSig(sigs, 0);
            ++sigsAdded;
        }

        if(!sigsInUse[4])
        {
            sigsInUse[4] = true;
            enableRefSig(sigs, 4);
            ++sigsAdded;
        }

        // force linked sigs to be included:

        // sig 7a with b,c and d
        if(sigsInUse[6] || sigsInUse[7] || sigsInUse[8] || sigsInUse[9])
        {
            for(int sig = 6; sig <= 9; ++sig)
            {
                if (sigsInUse[sig])
                    continue;

                sigsInUse[sig] = true;
                enableRefSig(sigs, sig);
                ++sigsAdded;
            }
        }

        // sig 10a with 10b
        if(sigsInUse[12] && !sigsInUse[13])
        {
            for(int sig = 12; sig <= 13; ++sig)
            {
                if (sigsInUse[sig])
                    continue;

                sigsInUse[sig] = true;
                enableRefSig(sigs, sig);
                ++sigsAdded;
            }
        }

        // sig 10a with 10b
        if(sigsInUse[20] && !sigsInUse[21])
        {
            for(int sig = 20; sig <= 21; ++sig)
            {
                if (sigsInUse[sig])
                    continue;

                sigsInUse[sig] = true;
                enableRefSig(sigs, sig);
                ++sigsAdded;
            }
        }

        return sigsAdded;
    }

    public static boolean fitCountsToRatios(
            int sampleId, final double[] counts, final double[] countsMargin, final List<double[]> ratiosCollection,
            double[] contribs, double minContribPerc)
    {
        // constraints are that the fitted counts cannot exceed the actual counts (plus the margin for noise)
        // if any single contrib drops below the required contribution percent, zero-out its sig

        int SCOL = 0; // the first and only column in the per-sample matrix, since only one set of data is handled at a time

        // validate inputs
        int bucketCount = counts.length;

        for(final double[] ratios : ratiosCollection)
        {
            if(ratios.length != bucketCount)
                return false;
        }

        int sigCount = ratiosCollection.size();

        double totalCount = sumVector(counts);

        // perform standard adjustments
        NmfMatrix w = new NmfMatrix(bucketCount, sigCount);
        NmfMatrix h = new NmfMatrix(sigCount, 1);

        for(int sig = 0; sig < sigCount; ++sig)
        {
            w.setCol(sig, ratiosCollection.get(sig));
            h.set(sig, SCOL, contribs[sig]);
        }

        NmfMatrix actualCounts = new NmfMatrix(bucketCount, 1);
        actualCounts.setCol(0, counts);

        NmfMatrix fittedCounts = new NmfMatrix(bucketCount, 1);

        // test out initial contributions and residuals
        w.multiply(h, fittedCounts, true);

        int iterations = 0;
        int iterationsCap = 50;
        int iterationsIncrement = 50;
        int maxIterations = 250;

        minContribPerc *= 0.9; // give a buffer in case a sig drops just below the min to allow it to resurface

        double currentCost = 0;
        double residuals = 0;
        double residualsPerc = 0;
        double prevCost = 0;
        double costChange = 0;
        double prevCostChange = 0;
        double prevResiduals = residuals;

        final double[][] fitData = fittedCounts.getData();

        List<Integer> zeroedSigs = Lists.newArrayList();

        while(iterations < iterationsCap)
        {
            // calc residuals
            prevCost = currentCost;
            prevResiduals = residuals;
            currentCost = actualCounts.sumDiffSq(fittedCounts);
            residuals = calcAbsDiffs(counts, fittedCounts.getCol(SCOL));
            residualsPerc = residuals / totalCount;

            if(Double.isNaN(currentCost) || Double.isInfinite(currentCost) || currentCost > 1e50)
            {
                LOGGER.warn("iter({}): invalid cost value: nan={} infinite={} max={}, exiting",
                        iterations, Double.isNaN(currentCost), Double.isInfinite(currentCost), currentCost > 1e50);

                return false;
            }

            if(residualsPerc < 0.01)
                break;

            if(iterations > 0)
            {
                prevCostChange = costChange;
                costChange = (prevCost - currentCost) / prevCost;

                if(abs(costChange) < 0.0001)
                    break;
            }

            // fit again
            NmfMatrix wt = w.transpose();
            NmfMatrix hAdj = wt.multiply(actualCounts);
            NmfMatrix hd = wt.multiply(fittedCounts);

            hAdj.scalarDivide(hd, true);
            h.scalarMultiply(hAdj);

            w.multiply(h, fittedCounts, true);

            // zero-out any tiny contributions
            for(int s = 0; s < sigCount; ++s)
            {
                if(zeroedSigs.contains(s))
                    continue;;

                double sigPerc = h.get(s, SCOL) / totalCount;
                if(sigPerc < minContribPerc)
                {
                    zeroedSigs.add(s);
                    h.set(s, SCOL, 0);

                    for(int b = 0; b < bucketCount; ++b)
                    {
                        w.set(b, s, 0);
                    }
                }
            }

            ++iterations;

            if(iterations >= iterationsCap && iterations <= maxIterations && costChange > 0.001) // 100 out of 100K, keep going if still progressing
            {
                iterationsCap += iterationsIncrement;
            }
        }

        // before beginning the next set of adjustments, check for any bucket exceeding the permitted range
        /* Ensure fitted counts are within the permitted noise range:
            - check each bucket's fitted vs actual count
            - if the fitted count exceeds the actual, calculate the required reduction per contributing sig for that bucket
         */

        double[] contribAdj = new double[sigCount];
        double[][] contribData = h.getData();
        final double[][] sigsData = w.getData();

        double contribTotal = h.sum();

        for(int b = 0; b < bucketCount; ++b)
        {
            double excessCount = fitData[b][SCOL] - (counts[b] + countsMargin[b]);

            if (excessCount <= 0)
                continue;

            double[] sigInvCosts = new double[sigCount];
            double invCostTotal = 0;

            for(int s = 0; s < sigCount; ++s)
            {
                // how much this sig contributed to this bucket
                double sigBucketContrib = contribData[s][SCOL] * sigsData[b][s];

                if (sigBucketContrib == 0)
                {
                    continue;
                }
                // the cost basis per unit
                double sigCostBasis = excessCount * (1 / sigsData[b][s]);

                // reversed percent allocation
                double invCost = 1 / sigCostBasis;

                invCostTotal += invCost;
                sigInvCosts[s] = invCost;
            }

            for(int s = 0; s < sigCount; ++s)
            {
                double sigAttribPerc = sigInvCosts[s] / invCostTotal;

                if(sigAttribPerc == 0)
                    continue;

                double sigAdjust = excessCount * sigAttribPerc * (1 / sigsData[b][s]);

                contribAdj[s] = max(contribAdj[s], sigAdjust);

                if(sigAdjust > contribData[s][SCOL] * 1.01)
                {
                    LOGGER.error("excess contrib adjust");
                }
            }
        }

        // set the final contributions minus any adjustments
        for(int s = 0; s < sigCount; ++s)
        {
            contribData[s][SCOL] -= contribAdj[s];
            contribs[s] = max(contribData[s][SCOL], 0);
        }

        double contribAdjustments = sumVector(contribAdj);

        if(contribAdjustments > 0)
        {
            w.multiply(h, fittedCounts, true);
            residuals = calcAbsDiffs(counts, fittedCounts.getCol(SCOL));
            residualsPerc = residuals / totalCount;
        }

        LOGGER.debug(String.format("sample(%d) sigs(%d -> %d) total(%s) contrib(init=%s less adj=%s) new residuals(%s perc=%.4f) iters(%d) costChg(%.4f -> %.4f)",
                sampleId, sigCount, sigCount - zeroedSigs.size(),
                sizeToStr(totalCount), sizeToStr(contribTotal), sizeToStr(contribAdjustments),
                sizeToStr(residuals), residualsPerc, iterations, prevCostChange, costChange));

        return true;
    }


}
