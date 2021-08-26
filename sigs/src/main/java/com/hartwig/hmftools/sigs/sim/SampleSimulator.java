package com.hartwig.hmftools.sigs.sim;

import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.MatrixUtils.createMatrixFromListData;
import static com.hartwig.hmftools.common.utils.MatrixUtils.writeMatrixData;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sigs.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sigs.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.getNewFile;
import static com.hartwig.hmftools.common.sigs.DataUtils.getPoissonRandom;
import static com.hartwig.hmftools.common.sigs.DataUtils.getPoissonRandomLarge;
import static com.hartwig.hmftools.common.utils.Matrix.extractNonZeros;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SampleSimulator
{
    private static final Logger LOGGER = LogManager.getLogger(SampleSimulator.class);

    private SimConfig mConfig;
    private GenericDataCollection mDataCollection;
    private String mOutputDir;
    private String mOutputFileId;

    // file of simulated signature parameters
    private final List<SimSigFactors> mSigFactors;

    // reference signatures used to generate bucket counts
    private Matrix mInputSignatures;

    // bucket counts per simulated sample
    private Matrix mOutputMatrix;

    // cache the allocation of sigs to each simulated sample
    private Matrix mOutputContributions;

    // random number caches
    private final double[] mPoissonDecimals;
    private final int[] mPoissonInts; // used for signature count distribution around the mean
    private int mPoissonIndex;

    // measure residuals from rounding the counts from fractions to integers
    private double mGrossResiduals;
    private double mNetResiduals;
    private int mGrossCountNoise;
    private int mNetCountNoise;

    private Random mRandom;

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
        mOutputContributions = null;
        mPoissonDecimals = new double[POISSON_DIST_SIZE];
        mPoissonInts = new int[POISSON_DIST_SIZE];
        mPoissonIndex = 0;
        mNetResiduals = 0;
        mGrossResiduals = 0;
        mGrossCountNoise = 0;
        mNetCountNoise = 0;

    }

    public void initialise(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mConfig = new SimConfig(cmd);

        mDataCollection = GenericDataLoader.loadFile(mConfig.SignaturesFilename);

        mRandom = mConfig.SeedRandom ? new Random(123456) : new Random();

        mInputSignatures = createMatrixFromListData(mDataCollection.getData());

        loadSignatureFactors(mConfig.SigFactorsFilename);

        LOGGER.info("loaded {} signature parameter sets:", mSigFactors.size());

        for(final SimSigFactors sigFactors : mSigFactors)
        {
            if (sigFactors.SigId-1 >= mInputSignatures.Cols)
            {
               LOGGER.error("sigFactor({}: {}) outside loaded signatures", sigFactors.SigId, sigFactors.Name);
            }

            LOGGER.info("{}", sigFactors.toString());
        }

        generatePoissonDist();
    }

    private void loadSignatureFactors(final String filename) {

        if (filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // read field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty data CSV file({})", filename);
                return;
            }

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                if(items.length != SimSigFactors.SIG_FACTOR_CSV_ITEM_COUNT)
                {
                    LOGGER.error("invalid signature parameter set: {}", line);
                    continue;
                }

                SimSigFactors sigFactors = new SimSigFactors(items);
                mSigFactors.add(sigFactors);
            }

            LOGGER.debug("loaded {} signature parameter sets", mSigFactors.size());

        }
        catch (IOException exception)
        {

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

        // remove any samples with zero counts
        mOutputMatrix = extractNonZeros(mOutputMatrix);

        if(mOutputMatrix.Cols < mConfig.SampleCount)
        {
            LOGGER.info("actualSampleCount({}) vs config() due to zero-count simulated counts", mOutputMatrix.Cols, mConfig.SampleCount);
        }

        perfCounter.stop();

        logBucketStats();
        writeSampleCounts();
        writeSampleContributions();

        perfCounter.logStats();
    }

    private void generateSimSampleCounts()
    {
        int buckets = mInputSignatures.Rows;
        int sampleCount = mConfig.SampleCount;

        mOutputMatrix = new Matrix(buckets, sampleCount);

        mOutputContributions = new Matrix(mSigFactors.size(), sampleCount);
        double[][] cData = mOutputContributions.getData();

        int sigIndex = 0;
        for(final SimSigFactors sigFactors : mSigFactors)
        {
            int samplesExcluded = 0;

            for (int n = 0; n < mConfig.SampleCount; ++n)
            {
                // determine whether this sig should have any presence in this sample
                double samSigProb = mRandom.nextDouble();

                if (samSigProb > sigFactors.SampleProbability)
                {
                    samplesExcluded++;
                    continue;
                }

                double timeFactor = mRandom.nextDouble() * sigFactors.TimeFactor;

                int variantCount = calcSigCount(sigFactors.MedianCount, timeFactor, sigFactors.RateFactor);

                LOGGER.debug("sig({}) sample({}) variantCount({})", sigFactors.SigId, n, variantCount);

                cData[sigIndex][n] = variantCount;

                setBucketCounts(variantCount, n, sigFactors.SigId, mOutputMatrix);
            }

            double samplePerc = samplesExcluded/(double)mConfig.SampleCount;

            LOGGER.debug(String.format("sig(%d: %s) samplesExcluded(%d asPerc=%.3f)",
                    sigFactors.SigId, sigFactors.Name, samplesExcluded, samplePerc));

            ++sigIndex;
        }
    }

    private int calcSigCount(int medianCount, double timeFactor, double rateFactor)
    {
        double variantCount = 0;

        // log-normal approach
        double meanNorm = log(medianCount);
        // double logNormCount = exp(meanNorm + rateFactor*random);

        // slow log-normal approach
        int poissonInt = getNextRandomInt();
        double logNormCount = exp(meanNorm + rateFactor * poissonInt);

        // straight random approach
//        double maxCount = medianCount * rateFactor * 10; // temporary
//        double simpleVC = medianCount + (maxCount-medianCount) * random;

        variantCount = logNormCount;
        // variantCount = simpleVC;

        variantCount *= timeFactor;

        return (int) round(variantCount);
    }

    private void setBucketCounts(int variantCount, int sampleIndex, int sigIndex, Matrix sampleCounts)
    {
        int sigLookupIndex = sigIndex-1; // since the signatures definition matrix is now zero-based

        // divide the variant count amongst the applicable buckets as per their defined proportions
        double sigTotal = 0;

        final double[][] sigData = mInputSignatures.getData();
        double[][] scData = sampleCounts.getData();

        for(int i = 0; i < mInputSignatures.Rows; ++i)
        {
            sigTotal += sigData[i][sigLookupIndex];
        }

        if(sigTotal <= 0)
            return;

        for(int i = 0; i < mInputSignatures.Rows; ++i)
        {
            double bucketPercent = sigData[i][sigLookupIndex] / sigTotal;

            double bucketSigRaw = bucketPercent * variantCount;
            int bucketSigCount = (int) round(bucketSigRaw);

            double fraction = bucketSigRaw - bucketSigCount;
            mNetResiduals += fraction;
            mGrossResiduals += abs(fraction);

            if(!mConfig.ApplyNoise)
            {
                scData[i][sampleIndex] += bucketSigCount;
                continue;
            }

            // optionally apply noise around the bucket counts
            int bucketSigCountAdj = applyNoise(bucketSigCount);

            scData[i][sampleIndex] += bucketSigCountAdj; // each sig's contribution is added

            int noiseDiff = bucketSigCountAdj - bucketSigCount;
            mGrossCountNoise += abs(noiseDiff);
            mNetCountNoise += noiseDiff;

//            if(bucketSigCount > 50) {

                LOGGER.debug(String.format("bucket(%d) count(%d raw=%.2f) adj(%d diff=%d asPerc=%.3f)",
                        i, bucketSigCount, bucketSigRaw, bucketSigCountAdj, noiseDiff,
                        bucketSigCount > 0 ? noiseDiff / (double) bucketSigCount : 0));
//            }
        }
    }

    private int applyNoise(int bucketCount)
    {
        if(bucketCount <= 1)
            return bucketCount;
        else if(bucketCount <= 100)
            return getPoissonRandom(bucketCount, mRandom);
        else
            return getPoissonRandomLarge(bucketCount, mRandom);
    }

    private void logBucketStats()
    {
        // calculate counts, min, max, mean and median values per bucket
        final double[][] scData = mOutputMatrix.getData();

        double totalCount = mOutputMatrix.sum();
        int sampleCount = mOutputMatrix.Cols;
        int medianIndex = sampleCount / 2; // not averaged for even sample counts

        // by bucket
        for (int i = 0; i < mOutputMatrix.Rows; ++i)
        {
            int total = 0;
            int min = 0;
            int max = 0;
            double median = 0;

            final double[] bucketCounts = mOutputMatrix.getRow(i);

            List<Integer> sortedIndices = getSortedVectorIndices(bucketCounts, true);

            for (int j = 0; j < mOutputMatrix.Cols; ++j)
            {
                int sampleIndex = sortedIndices.get(j);
                int count = (int)scData[i][sampleIndex];
                total += count;

                if(min == 0 || (count > 0 && count < min))
                    min = count;

                if(count > max)
                    max = count;

                if(j == medianIndex)
                {
                    median = count;
                }
            }

            double mean = total/(double)sampleCount;
            LOGGER.debug(String.format("bucket(%d) min(%d) max(%d) total(%d) mean(%.0f) median(%.0f) percOfTotal(%.4f)",
                    i, min, max, total, mean, median, total/totalCount));
        }

        if(totalCount > 0)
        {
            LOGGER.debug(String.format("rounding residuals: gross(%.1f perc=%.4f) net(%.1f perc=%.4f) vs total(%.0f)",
                    mGrossResiduals, mGrossResiduals / totalCount, mNetResiduals, mNetResiduals / totalCount, totalCount));

            if(mConfig.ApplyNoise)
            {
                LOGGER.debug(String.format("counts noise: gross(%d perc=%.4f) net(%d perc=%.4f) vs total(%.0f)",
                        mGrossCountNoise, mGrossCountNoise / totalCount, mNetCountNoise, mNetCountNoise / totalCount, totalCount));
            }
        }
    }

    public void writeSampleCounts()
    {
        try
        {
            final String filename = mOutputFileId + "_sim_sc.csv";
            BufferedWriter writer = getNewFile(mOutputDir, filename);

            int i = 0;
            for(; i < mConfig.SampleCount-1; ++i)
            {
                writer.write(String.format("Sample_%04d,", i));
            }
            writer.write(String.format("Sample_%04d", i));

            writer.newLine();

            writeMatrixData(writer, mOutputMatrix, true);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void writeSampleContributions()
    {
        try
        {
            final String filename = mOutputFileId + "_sim_contributions.csv";
            BufferedWriter writer = getNewFile(mOutputDir, filename);

            int i = 0;
            for(; i < mConfig.SampleCount-1; ++i)
            {
                writer.write(String.format("Sample_%d,", i));
            }
            writer.write(String.format("Sample_%d", i));

            writer.newLine();

            writeMatrixData(writer, mOutputContributions, false);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    private void generatePoissonDist()
    {
        // generate doubles between 0-1 based on a poisson distribution
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        SIG_LOGGER.info("running signature simulation");

        Options options = new Options();
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        SimConfig.addCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        SampleSimulator sampleSimulator = new SampleSimulator();
        sampleSimulator.initialise(cmd);
        sampleSimulator.run();

        SIG_LOGGER.info("signature simulation complete");

    }
}
