package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.data_analyser.DataAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getNewFile;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.getPoissonRandom;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.writeMatrixData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.data_analyser.loaders.GenericDataLoader;
import com.hartwig.hmftools.data_analyser.types.GenericDataCollection;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleSimulator {

    private static final Logger LOGGER = LogManager.getLogger(SampleSimulator.class);

    private SimConfig mConfig;
    private GenericDataCollection mDataCollection;
    private String mOutputDir;
    private String mOutputFileId;

    private List<SimSigFactors> mSigFactors;
    private NmfMatrix mOutputMatrix;
    private NmfMatrix mSignatures;
    private double[] mPoissonDecimals;
    private int[] mPoissonInts;
    private int mPoissonIndex;

    Random mRandom;

    private int SIG_FACTOR_CSV_ITEM_COUNT = 6;
    private int POISSON_DIST_SIZE = 1000;
    private int POISSON_LAMBDA = 15; // value around which the distribution is centred

    public SampleSimulator()
    {
        mOutputDir = "";
        mOutputFileId = "";
        mDataCollection = null;
        mSigFactors = Lists.newArrayList();
        mConfig = null;
        mRandom = null;
        mOutputMatrix = null;
        mPoissonDecimals = null;
        mPoissonInts = null;
        mPoissonIndex = 0;
    }

    public void initialise(final CommandLine cmd)
    {
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mConfig = new SimConfig(cmd);

        mDataCollection = GenericDataLoader.loadFile(mConfig.SignaturesFilename);

        mRandom = mConfig.SeedRandom ? new Random(123456) : new Random();

        mSignatures = DataUtils.createMatrixFromListData(mDataCollection.getData());

        loadSignatureFactors(mConfig.SigFactorsFilename);

        LOGGER.info("loaded {} signature parameter sets:", mSigFactors.size());
        for(final SimSigFactors sigFactors : mSigFactors) {

            if (sigFactors.SigId >= mSignatures.Cols) {

                LOGGER.error("sigFactor({}: {}) outside loaded signatures", sigFactors.SigId, sigFactors.Name);
            }

            LOGGER.info("{}", sigFactors.toString());
        }

        generatePoissonDist();
    }

    private void loadSignatureFactors(final String filename) {

        if (filename.isEmpty()) {
            return;
        }

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if (line == null) {
                LOGGER.error("Empty data CSV file({})", filename);
                return;
            }

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                if(items.length != SIG_FACTOR_CSV_ITEM_COUNT)
                {
                    LOGGER.error("invalid signature parameter set: {}", line);
                    continue;
                }

                SimSigFactors sigFactors = new SimSigFactors(items);
                mSigFactors.add(sigFactors);
            }

            LOGGER.debug("loaded {} signature parameter sets", mSigFactors.size());

        } catch (IOException exception) {

            LOGGER.error("failed to read data file({})", filename);
            return;
        }
    }

    public void run()
    {
        PerformanceCounter perfCounter = new PerformanceCounter("SampleSims");
        perfCounter.start();

        // testSampleCounts();
        generateSimSampleCounts();

        perfCounter.stop();
        perfCounter.logStats();

        logBucketStats();
        writeSampleCounts();
    }

    private void generateSimSampleCounts()
    {
        int buckets = mSignatures.Rows;
        int samples = mConfig.SampleCount;

        mOutputMatrix = new NmfMatrix(buckets, samples);

        for(final SimSigFactors sigFactors : mSigFactors) {
            int samplesExcluded = 0;

            for (int n = 0; n < mConfig.SampleCount; ++n) {
                // determine whether this sig should have any presence in this sample
                double samSigProb = mRandom.nextDouble();

                if (samSigProb > sigFactors.SampleProbability) {
                    samplesExcluded++;
                    continue;
                }

                double timeFactor = mRandom.nextDouble() * sigFactors.TimeFactor;

                int variantCount = calcSigCount(sigFactors.MedianCount, timeFactor, sigFactors.RateFactor);

                LOGGER.debug("sig({}) sample({}) variantCount({})", sigFactors.SigId, n, variantCount);

                setBucketCounts(variantCount, n, sigFactors.SigId, mOutputMatrix);
            }

            double samplePerc = samplesExcluded/(double)mConfig.SampleCount;

            LOGGER.debug(String.format("sig(%d: %s) samplesExcluded(%d asPerc=%.3f)",
                    sigFactors.SigId, sigFactors.Name, samplesExcluded, samplePerc));
        }
    }

    private int calcSigCount(int medianCount, double timeFactor, double rateFactor)
    {
        double variantCount = 0;
        // double random = getNextRandom();
        // double random = mRandom.nextDouble();

        // log-normal approach
        double meanNorm = log(medianCount);
        // double logNormCount = exp(meanNorm + rateFactor*random);

        // slow log-normal approach
        int poissonInt = getNextRandomInt();
        double logNormCount = exp(meanNorm + rateFactor*poissonInt);

//        int randomInt = mRandom.nextInt(20) - 10;
//        double logNormCount = exp(meanNorm + rateFactor*randomInt);

        // straight random approach
//        double maxCount = medianCount * rateFactor * 10; // temporary
//        double simpleVC = medianCount + (maxCount-medianCount) * random;

        variantCount = logNormCount;
        // variantCount = simpleVC;

        variantCount *= timeFactor;

        return (int) round(variantCount);
    }

    private void setBucketCounts(int variantCount, int sampleIndex, int sigIndex, NmfMatrix sampleCounts)
    {
        int sigLookupIndex = sigIndex-1; // since the signatures definititon matrix is now zero-based

        // divide the variant count amongst the applicable buckets as per their defined proportions
        double sigTotal = 0;

        final double[][] sigData = mSignatures.getData();
        double[][] scData = sampleCounts.getData();

        for(int i = 0; i < mSignatures.Rows; ++i)
        {
            sigTotal += sigData[i][sigLookupIndex];
        }

        if(sigTotal <= 0)
            return;

        for(int i = 0; i < mSignatures.Rows; ++i)
        {
            double bucketPercent = sigData[i][sigLookupIndex] / sigTotal;
            int bucketSigCount = (int) round(bucketPercent * variantCount);
            scData[i][sampleIndex] += bucketSigCount; // each sig's contribution is added
        }

        // optionally apply noise around the bucket counts
    }

    private void logBucketStats() {
        // calculate counts, min, max, mean and median values per signature and bucket
        final double[][] scData = mOutputMatrix.getData();

        // by bucket
        for (int i = 0; i < mOutputMatrix.Rows; ++i) {

            int total = 0;
            int min = 0;
            int max = 0;

            for (int j = 0; j < mOutputMatrix.Cols; ++j) {

                int count = (int)scData[i][j];
                total += count;

                if(min == 0 || (count > 0 && count < min))
                    min = count;

                if(count > max)
                    max = count;
            }

            double avg = total/(double)mOutputMatrix.Cols;
            LOGGER.debug(String.format("bucket(%d) min(%d) max(%d) total(%d) avg(%.0f)", i, min, max, total, avg));
        }
    }

    public void writeSampleCounts()
    {
        try
        {
            final String filename = mOutputFileId + "_sim_sample_counts.csv";
            BufferedWriter writer = getNewFile(mOutputDir, filename);

            int i = 0;
            for(; i < mConfig.SampleCount-1; ++i)
            {
                writer.write(String.format("Sample_%d,", i));
            }
            writer.write(String.format("Sample_%d", i));

            writer.newLine();

            writeMatrixData(writer, mOutputMatrix, true);

            writer.close();
        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    private void generatePoissonDist()
    {
        // generate doubles between 0-1 based on a poisson distributino
        mPoissonDecimals = new double[POISSON_DIST_SIZE];
        mPoissonInts = new int[POISSON_DIST_SIZE];

        int lambda = POISSON_LAMBDA;

        int[] freqCounts = new int[lambda*2+1];

        for(int i = 0; i < POISSON_DIST_SIZE; ++i)
        {
            int value = min(getPoissonRandom(lambda, mRandom),100);
            mPoissonInts[i] = value - POISSON_LAMBDA;
            mPoissonDecimals[i] = value/100.0;

            if(value < freqCounts.length)
                freqCounts[value] += 1;
        }
    }

    private double getNextRandom()
    {
        if(mPoissonIndex >= mPoissonDecimals.length)
            mPoissonIndex = 0;

        return mPoissonDecimals[mPoissonIndex++];
    }

    private int getNextRandomInt()
    {
        if(mPoissonIndex >= mPoissonInts.length)
            mPoissonIndex = 0;

        return mPoissonInts[mPoissonIndex++];
    }

    public void runTests()
    {
        testSampleCounts();
    }

    private void testSampleCounts()
    {
        int maxIterations = 100000;
        int lambda = 100;

        int[] freqCounts = new int[lambda*3];
//
//        maxIterations = matrix.Rows * matrix.Cols;

        int outOfBounds = 0;
        for(int i = 0; i < maxIterations; ++i)
        {
            int value = getPoissonRandom(lambda, mRandom);
            // LOGGER.debug(String.format("%d: value(%d)", i, value));

            if(value >= freqCounts.length)
                ++outOfBounds;
            else
                freqCounts[value] += 1;
        }

        for(int i = 0; i < freqCounts.length; ++i)
        {
            LOGGER.debug(String.format("POS: %d,%d", i, freqCounts[i]));
        }
    }

}
