package com.hartwig.hmftools.sigs.nmf;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.common.SigReporter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// fits a sample's counts to a pre-defined signature
// apply PCAWG fitting rules, which put constraints on how signatures are combined,
// and the expected mutational load for them

public class NmfSampleFitter
{
    private NmfConfig mConfig;
    private Matrix mRefSignatures;
    private Matrix mReferenceContribs;
    private Matrix mSampleCounts;
    private Matrix mAllContributions;

    private int[] mSigCountFrequency;

    private boolean mIsValid;

    // PCAWG rule constraints
    private static double SIG_EXCLUSION_CSS_DIFF = 0.01;
    private static double SIG_RE_INCLUSION_CSS_DIFF = 0.05;
    private static double MINOR_SIG_CONTRIBUTION_PERC = 0.001;

    private static final Logger LOGGER = LogManager.getLogger(NmfSampleFitter.class);

    public NmfSampleFitter(final NmfConfig config, final Matrix sampleCounts, final Matrix refSigs)
    {
        mConfig = config;
        mSampleCounts = sampleCounts;
        mAllContributions = new Matrix(refSigs.Cols, sampleCounts.Cols);
        mSigCountFrequency = new int[refSigs.Cols];
        mRefSignatures = refSigs;
        mReferenceContribs = null;
        mIsValid = true;
    }

    public void setRefContributions(final Matrix refContribs) { mReferenceContribs = refContribs; }
    public final Matrix getContributions() { return mAllContributions; }
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

    private boolean canSampleUseSig(int sampleId, int sig)
    {
        if(mReferenceContribs == null)
            return true;

        return mReferenceContribs.get(sig, sampleId) > 0;
    }

    private boolean fitSample(int sampleId)
    {
        int bucketCount = mSampleCounts.Rows;

        // prepare a matrix with only this sample's counts
        Matrix sampleMatrix = new Matrix(bucketCount, 1);

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
        Matrix reducedSigs = new Matrix(mRefSignatures);
        int currentSigCount = 0;

        // start with all in use unless below required ML threshold
        // final List<double[]> sigsCollection = Lists.newArrayList();
        for(int i = 0; i < refSigCount; ++i)
        {
            if(belowRequiredMutLoad(i, sampleCount) || !canSampleUseSig(sampleId, i))
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

            double newFitVsActualCss = calcCosineSim(newFit, sampleCounts);

            if(iterations > 0)
            {
                // assess the new contributions vs the actual counts, and exit if the change is greater than 0.01
                double cssDiff = lastFitVsActualCss - newFitVsActualCss;

                if (cssDiff > SIG_EXCLUSION_CSS_DIFF) // also apply this constraint when not applying other PCAWG fitting rules
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
                    if (leastIndex == -1 || newContribs[i] < leastContrib)
                    {
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

        double[] currentFit = nmfCalc.getFit().getCol(0);
        lastFitVsActualCss = calcCosineSim(currentFit, sampleCounts);

        LOGGER.debug(String.format("sample(%d) final run with %d sigs: css(%.4f) residuals(%.0f perc=%.3f) sampleCount(%.0f)",
                sampleId, currentSigCount, lastFitVsActualCss, nmfCalc.getTotalResiduals(),
                nmfCalc.getTotalResiduals()/sampleCount, sampleCount));

        if(mConfig.ApplyPcawgRules)
        {
            // finally try adding back in any signature which increases the CSS by the specified amount
            for (int i = 0; i < refSigCount; ++i)
            {
                if (sigsInUse[i] || !canSampleUseSig(sampleId, i))
                    continue;

                sigsInUse[i] = true;
                enableRefSig(reducedSigs, i);

                nmfCalc.setSignatures(reducedSigs);
                nmfCalc.performRun(sampleId);

                if (!nmfCalc.isValid())
                {
                    LOGGER.warn("NMF fit failed for sample({})", sampleId);
                    break;
                }

                final double[] newFit = nmfCalc.getFit().getCol(0);
                double newFitVsActualCss = calcCosineSim(newFit, sampleCounts);

                double cssDiff = newFitVsActualCss - lastFitVsActualCss;

                if (cssDiff >= SIG_RE_INCLUSION_CSS_DIFF)
                {
                    prevResiduals = nmfCalc.getTotalResiduals();
                    prevContribs = nmfCalc.getContributions().getCol(0);
                    ++currentSigCount;

                    LOGGER.debug(String.format("sample(%d) fitVsActualsCss(%.4f -> %.4f diff=%.4f) improved with new sig(%d) sigCount(%d)",
                            sampleId, lastFitVsActualCss, newFitVsActualCss, cssDiff, i, currentSigCount));
                }
                else
                {
                    reduceSigData(sigsInUse, reducedSigs, i);
                }

                lastFitVsActualCss = newFitVsActualCss;
            }
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

    private void reduceSigData(boolean[] sigsInUse, Matrix sigs, int sigIndex)
    {
        // disable a signature and set all contributions for it to zero
        sigsInUse[sigIndex] = false;

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex, 0);
        }
    }

    private void enableRefSig(Matrix sigs, int sigIndex)
    {
        final double[][] refData = mRefSignatures.getData();

        for (int i = 0; i < sigs.Rows; ++i)
        {
            sigs.set(i, sigIndex,  refData[i][sigIndex]);
        }
    }

    // PCAWG sig references
    // 0,   1,   2,   3,   4,   5,   6,    7,    8,    9,    10,  11,  12,    13,    14,   15,   16,   17,   18,   19,   20,    21,
    // SBS1,SBS2,SBS3,SBS4,SBS5,SBS6,SBS7a,SBS7b,SBS7c,SBS7d,SBS8,SBS9,SBS10a,SBS10b,SBS11,SBS12,SBS13,SBS14,SBS15,SBS16,SBS17a,SBS17b,
    // 22,   23,   24,   25,   26,   27,   28,    29,  30,   31,   32,    33,  34,   35,   36,   37,   38,   39,   40 etc
    // SBS18,SBS19,SBS20,SBS21,SBS22,SBS24,SBS26,SBS28,SBS30,SBS33,SBS34,SBS35,SBS36,SBS37,SBS38,SBS39,SBS40,SBSR1,SBSR2

    // PCAWG inclusion and exclusion rules
    public static final int PCAWG_SIG_1_AGE = 0;
    public static final int PCAWG_SIG_5_AGE = 4;
    public static final int PCAWG_SIG_2_AIDAPOBEC = 1;
    public static final int PCAWG_SIG_13_AIDAPOBEC = 16;
    public static final int PCAWG_SIG_7A_SKIN = 6;
    public static final int PCAWG_SIG_7B_SKIN = 7;
    public static final int PCAWG_SIG_7C_SKIN = 8;
    public static final int PCAWG_SIG_7D_SKIN = 9;
    public static final int PCAWG_SIG_10A = 12;
    public static final int PCAWG_SIG_10B = 13;
    public static final int PCAWG_SIG_17A = 20;
    public static final int PCAWG_SIG_17B = 21;

    private boolean belowRequiredMutLoad(int sig, double sampleCount)
    {
        if(!mConfig.ApplyPcawgRules)
            return false;

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

    private int addRequiredSigs(boolean[] sigsInUse, Matrix sigs)
    {
        if(!mConfig.ApplyPcawgRules)
            return 0;

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

        // 2 AID-APOBEC sigs (2 and 13)
        if(enforceSigPairInUse(PCAWG_SIG_2_AIDAPOBEC, PCAWG_SIG_13_AIDAPOBEC, sigsInUse, sigs))
        {
            ++sigsAdded;
        }

        // sig 10a with 10b
        if(enforceSigPairInUse(PCAWG_SIG_10A, PCAWG_SIG_10B, sigsInUse, sigs))
        {
            ++sigsAdded;
        }

        // sig 17a with 17b
        if(enforceSigPairInUse(PCAWG_SIG_17A, PCAWG_SIG_17B, sigsInUse, sigs))
        {
            ++sigsAdded;
        }

        return sigsAdded;
    }

    private boolean enforceSigPairInUse(int sig1, int sig2, boolean[] sigsInUse, Matrix sigs)
    {
        if(sigsInUse[sig1] == sigsInUse[sig2])
            return false;

        if (sigsInUse[sig1])
        {
            sigsInUse[sig2] = true;
            enableRefSig(sigs, sig2);
        }
        else
        {
            sigsInUse[sig1] = true;
            enableRefSig(sigs, sig1);
        }

        return true;
    }

}
