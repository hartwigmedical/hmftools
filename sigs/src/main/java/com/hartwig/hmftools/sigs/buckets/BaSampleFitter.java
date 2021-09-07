package com.hartwig.hmftools.sigs.buckets;

import static com.hartwig.hmftools.common.sigs.NoiseCalcs.calcRangeValue;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.BA_MAX_NOISE_ALLOC_PERCENT;
import static com.hartwig.hmftools.sigs.buckets.BaConfig.MAX_NOISE_ALLOC_PERCENT;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BaSampleFitter
{
    private final Matrix mSignatures;
    private final Matrix mContrbutions;
    private final Matrix mSampleCounts;
    private final List<String> mSampleIds;

    private final double mMinSigPercent;
    private final double mNoiseProbability;

    private final SampleSigContribOptimiser mSigContribOptimiser;

    private static final double DEFAULT_MIN_SIG_PERCENT = 0.01;
    private static final String MIN_SIG_PERCENT = "min_sig_percent";

    public static final double DEFAULT_NOISE_PROB = 1e-4;
    public static final String NOISE_PROB = "noise_prob";

    public static final int DEFAULT_MIN_NOISE_COUNT = 5;

    private static final Logger LOGGER = LogManager.getLogger(BaSampleFitter.class);

    public BaSampleFitter(final Matrix sampleCounts, final List<String> sampleIds, final Matrix signatures, final CommandLine cmd)
    {
        mSignatures = signatures;
        mSignatures.cacheTranspose();

        mSampleCounts = sampleCounts;
        mSampleCounts.cacheTranspose();
        mSampleIds = sampleIds;

        mContrbutions = new Matrix(signatures.Cols, sampleCounts.Cols);

        mMinSigPercent = Double.parseDouble(cmd.getOptionValue(MIN_SIG_PERCENT, String.valueOf(DEFAULT_MIN_SIG_PERCENT)));
        mNoiseProbability = Double.parseDouble(cmd.getOptionValue(NOISE_PROB, String.valueOf(DEFAULT_NOISE_PROB)));

        // set the max of a sample's total which can be allocated to noise, effecively setting a cap on
        MAX_NOISE_ALLOC_PERCENT = Double.parseDouble(cmd.getOptionValue(BA_MAX_NOISE_ALLOC_PERCENT, String.valueOf(MAX_NOISE_ALLOC_PERCENT)));

        mSigContribOptimiser = new SampleSigContribOptimiser(signatures.Rows, false, 1.0);
    }

    public final Matrix getContributions() { return mContrbutions; }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(MIN_SIG_PERCENT, true, "Min percent of sample's total count to allocate to signature");
        options.addOption(NOISE_PROB, true, "Poisson probability for noise per buckt");
        options.addOption(BA_MAX_NOISE_ALLOC_PERCENT, true, "Max percent of sample's total count to allocate to noise");
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

        final List<double[]> ratiosCollection = Lists.newArrayList();
        final List<Integer> sigIds = Lists.newArrayList();

        final Map<Integer,Integer> rangeCache = Maps.newHashMap();

        for(int s = 0; s < mSampleCounts.Cols; ++s)
        {
            SampleData sample = new SampleData(s);

            if(mSampleIds.size() == mSampleCounts.Cols)
            {
                sample.setSampleName(mSampleIds.get(s));
            }

            final double[] sampleCounts = mSampleCounts.getCol(s);
            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal == 0)
                continue;

            sample.setBucketCounts(sampleCounts);

            for(int b = 0; b < bucketCount; ++b)
            {
                if(sampleCounts[b] > 0)
                {
                    noiseCounts[b] = calcRangeValue(rangeCache, (int)sampleCounts[b], mNoiseProbability, DEFAULT_MIN_NOISE_COUNT, true);
                }
                else
                {
                    noiseCounts[b] = DEFAULT_MIN_NOISE_COUNT;
                }
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

            mSigContribOptimiser.initialise(sample, ratiosCollection, mMinSigPercent, 1);
            mSigContribOptimiser.setSigIds(sigIds);
            // mSigContribOptimiser.setLogVerbose(mConfig.logSample(sample.Id));

            // LOGGER.debug("sample({}) fitting with {} sigs", sample.Id, sigIds.size());

            boolean validCalc = mSigContribOptimiser.fitToSample();

            if (!validCalc)
            {
                LOGGER.error("sample({}) sig fit failed", sample.Id);
                return;
            }

            // if all ok, allocate each contribution to the sample
            double[] sigContribs = mSigContribOptimiser.getContribs();

            for(int j = 0; j < sigIds.size(); ++j)
            {
                int sigIndex = sigIds.get(j);
                mContrbutions.set(sigIndex, s, sigContribs[j]);
            }
        }
    }
}
