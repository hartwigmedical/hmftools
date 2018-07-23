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

    private boolean mIsValid;

    private static double SIG_CONTRIBUTION_SIMILARITY = 0.99;


    public NmfSampleFitter(final NmfConfig config, final NmfMatrix sampleCounts, final NmfMatrix refSigs)
    {
        mConfig = config;
        mSampleCounts = sampleCounts;
        mAllContributions = new NmfMatrix(refSigs.Cols, sampleCounts.Cols);
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

        // start with all in use
        for(int i = 0; i < refSigCount; ++i)
        {
            sigsInUse[i] = true;
        }

        NmfCalculator nmfCalc = new NmfCalculator(sampleCounts, mConfig);

        nmfCalc.setSigCount(refSigCount);

        double[] prevContribs = new double[refSigCount];
        int currentSigCount = refSigCount;
        double prevResiduals = 0;
        int iterations = 0;

        while(currentSigCount > 1)
        {
            nmfCalc.setSignatures(reducedSigs);
            nmfCalc.performRun(sampleId);

            if(!nmfCalc.isValid())
            {
                LOGGER.warn("NMF fit failed for sample({})", sampleId);
                return false;
            }

            final double[] newContribs = nmfCalc.getContributions().getCol(0);

            if(iterations > 0) {

                // assess the new contributions vs the prevous set, and exit if they're no longer similar
                double sigCss = calcCSS(newContribs, prevContribs);

                if (sigCss < SIG_CONTRIBUTION_SIMILARITY) {

                    LOGGER.debug(String.format("sample(%d) sig contrib css(%.4f) with sigCount(%d), exiting",
                            sampleId, sigCss, currentSigCount));
                    break;
                }
            }

            copyVector(newContribs, prevContribs);
            prevResiduals = nmfCalc.getTotalResiduals();

            // remove the least contributing signature
            int leastIndex = 0;
            double leastContrib = 0;

            for (int i = 0; i < refSigCount; ++i) {

                if(!sigsInUse[i])
                    continue;

                if (i == 0 || prevContribs[i] < leastContrib) {
                    leastIndex = i;
                    leastContrib = prevContribs[i];
                }
            }

            LOGGER.debug(String.format("sample(%d) remove least sig(%d contrib=%.4g)", sampleId, leastIndex, leastContrib));
            sigsInUse[leastIndex] = false;

            // now set that sig to zeros so it plays no part in the next fit
            for (int i = 0; i < reducedSigs.Rows; ++i) {
                reducedSigs.set(i, leastIndex, 0); // negligible but for NMF cannot be zero
            }

            --currentSigCount;

            ++iterations;
        }

        // use the previous fit's contributions ie when still had a sufficiently high CSS
        double[][] allContribData = mAllContributions.getData();

        for(int i = 0; i < refSigCount; ++i)
        {
            allContribData[i][sampleId] = prevContribs[i];
        }

        LOGGER.debug(String.format("sample(%d) residualsPerc(%.4f) sigCount(%d)",
                sampleId, prevResiduals / sampleCount, currentSigCount));

        return true;
    }

}
