package com.hartwig.hmftools.sig_analyser.buckets;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.SigMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BaSampleFitter
{
    private final SigMatrix mSignatures;
    private final SigMatrix mContrbutions;
    private final SigMatrix mSampleCounts;

    private final double mMinSigPercent;
    private static double DEFAULT_MIN_SIG_PERCENT = 0.01;
    private static final String MIN_SIG_PERCENT = "min_sig_percent";

    private final double mNoisePercent;
    private static double DEFAULT_NOISE_PERCENT = 0.05;
    private static final String NOISE_PERCENT = "noise_percent";
    private static int MIN_ABSOLUTE_NOISE = 5;

    private static final Logger LOGGER = LogManager.getLogger(BaSampleFitter.class);

    public BaSampleFitter(final SigMatrix sampleCounts, final SigMatrix signatures, final CommandLine cmd)
    {
        mSignatures = signatures;
        mSignatures.cacheTranspose();

        mSampleCounts = sampleCounts;
        mSampleCounts.cacheTranspose();

        mContrbutions = new SigMatrix(signatures.Cols, sampleCounts.Cols);

        mMinSigPercent = Double.parseDouble(cmd.getOptionValue(MIN_SIG_PERCENT, String.valueOf(DEFAULT_MIN_SIG_PERCENT)));
        mNoisePercent = Double.parseDouble(cmd.getOptionValue(NOISE_PERCENT, String.valueOf(DEFAULT_NOISE_PERCENT)));
    }

    public final SigMatrix getContributions() { return mContrbutions; }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(MIN_SIG_PERCENT, true, "Min percent of sample's total count to allocate to signature");
        options.addOption(NOISE_PERCENT, true, "Min percent of sample's total count to allocate to signature");
    }

    public void fitAllSamples()
    {
        // test each signature against the sample and then perform an fit optimisation with any which could contribute a sufficient amount
        int bucketCount = mSampleCounts.Rows;
        int sigCount = mSignatures.Cols;

        double[] noiseCounts = new double[bucketCount];
        double[] allocCounts = new double[bucketCount];
        double[] emptyBucketData = new double[bucketCount];

        final List<Integer> allBuckets = Lists.newArrayList();

        for(int b = 0; b < bucketCount; ++b)
        {
            allBuckets.add(b);
        }

        SigContribOptimiser sigContribOptimiser = new SigContribOptimiser(bucketCount, false, 1.0);
        List<double[]> ratiosCollection = Lists.newArrayList();
        List<Integer> sigIds = Lists.newArrayList();

        for(int s = 0; s < mSampleCounts.Cols; ++s)
        {
            SampleData sample = new SampleData(s);

            final double[] sampleCounts = mSampleCounts.getCol(s);
            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal == 0)
                continue;

            sample.setBucketCounts(sampleCounts);

            double minNoise = Arrays.stream(sampleCounts).filter(x -> x > 0).min().orElse(0);
            minNoise = minNoise > 0 ? min(minNoise, MIN_ABSOLUTE_NOISE) : MIN_ABSOLUTE_NOISE;

            for(int b = 0; b < bucketCount; ++b)
            {
                if(sampleCounts[b] > 0)
                    noiseCounts[b] = floor(mNoisePercent * sampleCounts[b]);
                else
                    noiseCounts[b] = minNoise;
            }

            sample.setElevatedBucketCounts(sample.getBucketCounts(), noiseCounts);

            ratiosCollection.clear();
            sigIds.clear();

            for(int i = 0; i < sigCount; ++i)
            {
                double[] sigRatios = mSignatures.getCol(i);
                double allocTotal = sample.getPotentialUnallocCounts(sigRatios, allBuckets, emptyBucketData, allocCounts);

                if (allocTotal / sampleTotal < mMinSigPercent)
                    continue;

                // add the sig's data to the optimiser
                ratiosCollection.add(sigRatios);
                sigIds.add(i);
            }

            sigContribOptimiser.initialise(sample, ratiosCollection, mMinSigPercent, 1);
            sigContribOptimiser.setSigIds(sigIds);
            // sigContribOptimiser.setLogVerbose(mConfig.logSample(sample.Id));

            // LOGGER.debug("sample({}) fitting with {} sigs", sample.Id, sigIds.size());

            boolean validCalc = sigContribOptimiser.fitToSample();

            if (!validCalc)
            {
                LOGGER.error("sample({}) sig fit failed", sample.Id);
                return;
            }

            // if all ok, allocate each contribution to the sample
            double[] sigContribs = sigContribOptimiser.getContribs();

            for(int j = 0; j < sigIds.size(); ++j)
            {
                int sigIndex = sigIds.get(j);
                mContrbutions.set(sigIndex, s, sigContribs[j]);
            }
        }
    }
}
