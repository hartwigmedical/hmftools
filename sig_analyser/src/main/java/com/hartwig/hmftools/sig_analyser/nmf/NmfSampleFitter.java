package com.hartwig.hmftools.sig_analyser.nmf;


import static com.hartwig.hmftools.sig_analyser.common.CosineSim.calcCSS;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.copyVector;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sig_analyser.buckets.SigContribOptimiser;
import com.hartwig.hmftools.sig_analyser.common.SigReporter;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// fits a sample's counts to a pre-defined signature
// apply PCAWG fitting rules, which put constraints on how signatures are combined,
// and the expected mutational load for them

public class NmfSampleFitter
{
   private NmfConfig mConfig;
    private SigMatrix mRefSignatures;
    private SigMatrix mSampleCounts;
    private SigMatrix mAllContributions;

    private int[] mSigCountFrequency;

    private boolean mIsValid;

    private static double MAX_FIT_CSS_DIFF = 0.01;
    private static double MINOR_SIG_CONTRIBUTION_PERC = 0.001;

    private static final Logger LOGGER = LogManager.getLogger(NmfSampleFitter.class);

    public NmfSampleFitter(final NmfConfig config, final SigMatrix sampleCounts, final SigMatrix refSigs)
    {
        mConfig = config;
        mSampleCounts = sampleCounts;
        mAllContributions = new SigMatrix(refSigs.Cols, sampleCounts.Cols);
        mSigCountFrequency = new int[refSigs.Cols];
        mRefSignatures = refSigs;
        mIsValid = true;
    }

    public final SigMatrix getContributions() { return mAllContributions; }
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

        // prepare a matrix with only this sample's counts
        SigMatrix sampleMatrix = new SigMatrix(bucketCount, 1);

        final double[] sampleCounts = mSampleCounts.getCol(sampleId);
        // final double[] sampleNoise = new double[bucketCount];
        double[][] sData = sampleMatrix.getData();
        double sampleCount = 0;

        for(int i = 0; i < bucketCount; ++i)
        {
            sData[i][0] = sampleCounts[i];
            sampleCount += sampleCounts[i];
        }

        int refSigCount = mRefSignatures.Cols;

        // keep track of which sigs are active
        boolean[] sigsInUse = new boolean[refSigCount];

        // inactive sigs will be zeroed out in the sigs matrix given to the NMF calculator
        SigMatrix reducedSigs = new SigMatrix(mRefSignatures);
        int currentSigCount = 0;

        // start with all in use unless below required ML threshold
        // final List<double[]> sigsCollection = Lists.newArrayList();
        for(int i = 0; i < refSigCount; ++i)
        {
            if(belowRequiredMutLoad(i, sampleCount))
            {
                reduceSigData(sigsInUse, reducedSigs, i);
            }
            else
            {
                sigsInUse[i] = true;
                // sigsCollection.add(reducedSigs.getCol(i));
                ++currentSigCount;
            }
        }

        NmfCalculator nmfCalc = new NmfCalculator(sampleMatrix, mConfig);

        nmfCalc.setSigCount(refSigCount);

        // SigContribOptimiser sigOptim = new SigContribOptimiser(bucketCount, false, 1.0);
        // sigOptim.initialise(sampleId, sampleCounts, sampleNoise, sigsCollection, 0.001, 0);

        double[] prevContribs = new double[refSigCount];
        double prevResiduals = 0;
        double lastFitVsActualCss = 0;
        int lastSigRemoved = 0;
        int iterations = 0;

        while(currentSigCount >= 1)
        {
            nmfCalc.setSignatures(reducedSigs);
            nmfCalc.performRun(sampleId);

            // boolean calcOk = sigOptim.fitToSample();

            if(!nmfCalc.isValid())
            {
                LOGGER.warn("NMF fit failed for sample({})", sampleId);
                return false;
            }

            final double[] newContribs = nmfCalc.getContributions().getCol(0);
            final double[] newFit = nmfCalc.getFit().getCol(0);

            double newFitVsActualCss = calcCSS(newFit, sampleCounts);

            if(iterations > 0)
            {
                // assess the new contributions vs the actual counts, and exit if the change is greater than 0.01
                double cssDiff = lastFitVsActualCss - newFitVsActualCss;

                if (cssDiff > MAX_FIT_CSS_DIFF)
                {
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

            for (int i = 0; i < refSigCount; ++i)
            {
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

                if(minorSigsRemoved == 0)
                {
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

    private void reduceSigData(boolean[] sigsInUse, SigMatrix sigs, int sigIndex)
    {
        // disable a signature and set all contributions for it to zero
        sigsInUse[sigIndex] = false;

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex, 0);
        }
    }

    private void enableRefSig(SigMatrix sigs, int sigIndex)
    {
        final double[][] refData = mRefSignatures.getData();

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex,  refData[i][sigIndex]);
        }
    }

    // PCAWG inclusion and exclusion rules
    public final static int PCAWG_SIG_1_AGE = 0;
    public final static int PCAWG_SIG_5_AGE = 4;
    public final static int PCAWG_SIG_7A_SKIN = 6;
    public final static int PCAWG_SIG_7B_SKIN = 7;
    public final static int PCAWG_SIG_7C_SKIN = 8;
    public final static int PCAWG_SIG_7D_SKIN = 9;
    public final static int PCAWG_SIG_10A = 12;
    public final static int PCAWG_SIG_10B = 13;
    public final static int PCAWG_SIG_17A = 20;
    public final static int PCAWG_SIG_17B = 21;

    private boolean belowRequiredMutLoad(int sig, double sampleCount)
    {
        // some signatures can only be applied for samples above a certain mutational load
        switch(sig)
        {
            case PCAWG_SIG_10A:
            case PCAWG_SIG_10B:
                return sampleCount < 1e5;

            case 5: // sig 6
            case 17: // sig 14
            case 18: // sig 15
            case 24: // 20
            case 25: // 21
            case 30: // 26
            case 48:
                return sampleCount < 1e4;

            default:
                return false;
        }

    }

    private int addRequiredSigs(boolean[] sigsInUse, SigMatrix sigs)
    {
        // force some signatures to be included: currently 1 and 5
        int sigsAdded = 0;

        if(!sigsInUse[PCAWG_SIG_1_AGE])
        {
            sigsInUse[PCAWG_SIG_1_AGE] = true;
            enableRefSig(sigs, PCAWG_SIG_1_AGE);
            ++sigsAdded;
        }

        if(!sigsInUse[PCAWG_SIG_5_AGE])
        {
            sigsInUse[PCAWG_SIG_5_AGE] = true;
            enableRefSig(sigs, PCAWG_SIG_5_AGE);
            ++sigsAdded;
        }

        // force linked sigs to be included:

        // sig 7a with b,c and d
        if(sigsInUse[PCAWG_SIG_7A_SKIN] || sigsInUse[PCAWG_SIG_7B_SKIN] || sigsInUse[PCAWG_SIG_7C_SKIN] || sigsInUse[PCAWG_SIG_7D_SKIN])
        {
            for(int sig = PCAWG_SIG_7A_SKIN; sig <= PCAWG_SIG_7D_SKIN; ++sig)
            {
                if (sigsInUse[sig])
                    continue;

                sigsInUse[sig] = true;
                enableRefSig(sigs, sig);
                ++sigsAdded;
            }
        }

        // sig 10a with 10b
        if(sigsInUse[PCAWG_SIG_10A] && !sigsInUse[PCAWG_SIG_10B])
        {
            for(int sig = PCAWG_SIG_10A; sig <= PCAWG_SIG_10B; ++sig)
            {
                if (sigsInUse[sig])
                    continue;

                sigsInUse[sig] = true;
                enableRefSig(sigs, sig);
                ++sigsAdded;
            }
        }

        // sig 17a with 17b
        if(sigsInUse[PCAWG_SIG_17A] && !sigsInUse[PCAWG_SIG_17B])
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

}
