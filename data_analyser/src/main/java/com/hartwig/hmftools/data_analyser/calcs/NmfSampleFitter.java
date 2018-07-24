package com.hartwig.hmftools.data_analyser.calcs;

import static com.hartwig.hmftools.data_analyser.calcs.CosineSim.calcCSS;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;

import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

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

    private static double SIG_CONTRIBUTION_SIMILARITY = 0.99;
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

            LOGGER.debug("assigned sig count frequency({}) = {}", i, mSigCountFrequency[i]);
        }

        SigReporter sigReporter = new SigReporter(mSampleCounts, mRefSignatures, mAllContributions, mRefSignatures, mRefSignatures, mConfig);
        sigReporter.runAnalysis();
    }

    private boolean fitSample(final int sampleId)
    {
        int bucketCount = mSampleCounts.Rows;

        NmfMatrix sampleCounts = new NmfMatrix(bucketCount, 1);

        final double[][] allData = mSampleCounts.getData();
        double[][] sData = sampleCounts.getData();
        double sampleCount = 0;

        for(int i = 0; i < bucketCount; ++i)
        {
            sData[i][0] = allData[i][sampleId];
            sampleCount += allData[i][sampleId];
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

        NmfCalculator nmfCalc = new NmfCalculator(sampleCounts, mConfig);

        nmfCalc.setSigCount(refSigCount);

        double[] prevContribs = new double[refSigCount];
        double prevResiduals = 0;
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
            double contribsVsLastCss = 0;

            if(iterations > 0) {

                // assess the new contributions vs the prevous set, and exit if they're no longer similar
                contribsVsLastCss = calcCSS(newContribs, prevContribs);

                if (contribsVsLastCss < SIG_CONTRIBUTION_SIMILARITY) {

                    LOGGER.debug(String.format("sample(%d) sig contrib css(%.4f) with sigCount(%d), exiting",
                            sampleId, contribsVsLastCss, currentSigCount));

                    // restore to last sig count prior to drop-off
                    ++currentSigCount;

                    break;
                }
            }

            copyVector(newContribs, prevContribs);
            prevResiduals = nmfCalc.getTotalResiduals();

            if(currentSigCount == 1)
            {
                LOGGER.debug(String.format("sample(%d) sig contrib css(%.4f) single sig, exiting",
                        sampleId, contribsVsLastCss));
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

                LOGGER.debug(String.format("sample(%d) remove least sig(%d contrib=%.1f perc=%.3f), cssVsLast(%.4f) remaining(%d)",
                        sampleId, leastIndex, leastContrib, leastContrib/sampleCount, contribsVsLastCss, currentSigCount));

                reduceSigData(sigsInUse, reducedSigs, leastIndex);
            }
            else
            {
                LOGGER.debug(String.format("sample(%d) no more sigs removed, count(%d)", sampleId, currentSigCount));
                break;
            }

            ++iterations;
        }

        // use the previous fit's contributions ie when still had a sufficiently high CSS
        double[][] allContribData = mAllContributions.getData();

        for(int i = 0; i < refSigCount; ++i)
        {
            allContribData[i][sampleId] = prevContribs[i];
        }

        LOGGER.debug(String.format("sample(%d) residualsPerc(%.0f of %.0f perc=%.4f) sigCount(%d)",
                sampleId, prevResiduals, sampleCount, prevResiduals/sampleCount, currentSigCount));

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

        for (int i = 0; i < sigs.Rows; ++i) {
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

}
