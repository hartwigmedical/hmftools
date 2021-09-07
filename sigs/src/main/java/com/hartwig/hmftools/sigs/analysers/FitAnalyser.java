package com.hartwig.hmftools.sigs.analysers;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.NoiseCalcs.calcRangeValue;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_UNALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.MatrixUtils.writeMatrixData;
import static com.hartwig.hmftools.common.utils.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sigs.buckets.BaSampleFitter.DEFAULT_MIN_NOISE_COUNT;
import static com.hartwig.hmftools.sigs.buckets.BaSampleFitter.DEFAULT_NOISE_PROB;
import static com.hartwig.hmftools.sigs.buckets.BaSampleFitter.NOISE_PROB;
import static com.hartwig.hmftools.sigs.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_IDS;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.ExpectationMaxFit;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.sigs.buckets.BaSampleFitter;
import com.hartwig.hmftools.sigs.common.CommonUtils;
import com.hartwig.hmftools.sigs.fitter.FitMethod;
import com.hartwig.hmftools.sigs.nmf.NmfConfig;
import com.hartwig.hmftools.sigs.nmf.NmfSampleFitter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class FitAnalyser
{
    // config
    private final List<FitMethod> mFitMethods;
    private final String mOutputDir;
    private final String mOutputId;
    private final CommandLine mCmdLineArgs;

    private final Matrix mSampleCounts;
    private final List<String> mSampleIds;
    private final Matrix mSignatures;
    private final List<String> mSigNames;

    private final Map<Integer,Integer> mNoiseRangeMap;
    private final double mNoiseProbability;

    private BufferedWriter mSampleResultsWriter;
    private BufferedWriter mResidualsWriter;

    // config strings
    private static final String FIT_METHODS = "fit_methods";
    private static final String SIGNATURES_FILE = "signatures_file";

    public FitAnalyser(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mCmdLineArgs = cmd;
        mFitMethods = Lists.newArrayList();

        final String[] fitMethods = cmd.getOptionValue(FIT_METHODS).split(";", -1);
        Arrays.stream(fitMethods).map(x -> FitMethod.valueOf(x)).forEach(x -> mFitMethods.add(x));

        mSampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_IDS))
        {
            mSampleIds.addAll(Arrays.stream(cmd.getOptionValue(SAMPLE_IDS).split(";", -1)).collect(Collectors.toList()));

            final List<String> sampleIds = Lists.newArrayList();
            Matrix allCounts = loadMatrixDataFile(cmd.getOptionValue(SAMPLE_COUNTS_FILE), sampleIds);

            mSampleCounts = new Matrix(allCounts.Rows, mSampleIds.size());

            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                for(int j = 0; j < sampleIds.size(); ++j)
                {
                    if(sampleIds.get(j).equals(mSampleIds.get(i)))
                        mSampleCounts.setCol(i, allCounts.getCol(j));
                }
            }
        }
        else
        {
            mSampleCounts = loadMatrixDataFile(cmd.getOptionValue(SAMPLE_COUNTS_FILE), mSampleIds);
        }

        mSampleCounts.cacheTranspose();

        mSigNames = Lists.newArrayList();
        mSignatures = loadMatrixDataFile(cmd.getOptionValue(SIGNATURES_FILE), mSigNames);
        mSignatures.cacheTranspose();

        mNoiseRangeMap = Maps.newHashMap();
        mNoiseProbability = Double.parseDouble(cmd.getOptionValue(NOISE_PROB, String.valueOf(DEFAULT_NOISE_PROB)));

        mSampleResultsWriter = null;
        mResidualsWriter = null;
        intialiseOutputFiles();
    }

    public void run()
    {
        if(mFitMethods.isEmpty() || mOutputDir == null)
            return;

        for(FitMethod fitMethod : mFitMethods)
        {
            SIG_LOGGER.info("method({}) fitting {} samples with {} signatures", fitMethod, mSampleCounts.Cols, mSignatures.Cols);

            Matrix sampleContribs = null;

            switch(fitMethod)
            {
                case LEAST_SQUARES:
                    sampleContribs = fitWithLeastSquares();
                    break;

                case EXPECTATIONS_MAX:
                    sampleContribs = fitWithExpectationsMax();
                    break;

                case SIG_OPTIMISER:
                    sampleContribs = fitWithSigOptimiser();
                    break;

                case NMF:
                    sampleContribs = fitWwithNMF();
                    break;

                default:
                    break;
            }

            processFitResults(fitMethod, sampleContribs);
            writeSigContributions(fitMethod, sampleContribs);
        }

        closeBufferedWriter(mSampleResultsWriter);
        closeBufferedWriter(mResidualsWriter);
    }

    private Matrix fitWithLeastSquares()
    {
        final Matrix sampleContribs = new Matrix(mSignatures.Cols, mSampleCounts.Cols);

        LeastSquaresFit lsqFit = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);

        for(int i = 0; i < mSampleCounts.Cols; ++i)
        {
            final double[] sampleCounts = mSampleCounts.getCol(i);
            lsqFit.initialise(mSignatures.getData(), sampleCounts);
            lsqFit.solve();

            final double[] sigAllocs = lsqFit.getContribs();
            sampleContribs.setCol(i, sigAllocs);
        }

        return sampleContribs;
    }

    private Matrix fitWithExpectationsMax()
    {
        final Matrix sampleContribs = new Matrix(mSignatures.Cols, mSampleCounts.Cols);

        for(int i = 0; i < mSampleCounts.Cols; ++i)
        {
            final double[] sampleCounts = mSampleCounts.getCol(i);
            final double[] sigAllocs = ExpectationMaxFit.performFit(sampleCounts, mSignatures, 0.001, 100);
            sampleContribs.setCol(i, sigAllocs);
        }

        return sampleContribs;
    }

    private Matrix fitWithSigOptimiser()
    {
        final Matrix sampleContribs = new Matrix(mSignatures.Cols, mSampleCounts.Cols);

        BaSampleFitter sampleFitter = new BaSampleFitter(mSampleCounts, mSampleIds, mSignatures, mCmdLineArgs);
        sampleFitter.fitAllSamples();
        sampleContribs.setData(sampleFitter.getContributions().getData());

        return sampleContribs;
    }

    private Matrix fitWwithNMF()
    {
        final Matrix sampleContribs = new Matrix(mSignatures.Cols, mSampleCounts.Cols);

        NmfConfig nmfConfig = new NmfConfig(1, 100);

        NmfSampleFitter nmfFitter = new NmfSampleFitter(nmfConfig, mSampleCounts, mSignatures);

        nmfFitter.fitSamples();

        sampleContribs.setData(nmfFitter.getContributions().getData());
        return sampleContribs;
    }

    private void writeSigContributions(final FitMethod fitMethod, final Matrix sampleContribs)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputDir + "FIT_" + fitMethod.toString() + "_contribs.csv", false);
            writeMatrixData(writer, mSampleIds, sampleContribs, false);
        }
        catch (final IOException e)
        {
            SIG_LOGGER.error("error writing sig contributions: {}", e.toString());
        }
    }

    private void processFitResults(final FitMethod fitMethod, final Matrix sampleContribs)
    {
        double totalCounts = 0;
        double totalAlloc = 0;
        double totalResiduals = 0;
        double totalResidualsExcess = 0;

        int sampleCount = mSampleCounts.Cols;

        for(int i = 0; i < sampleCount; ++i)
        {
            final String sampleId = mSampleIds.get(i);
            final double[] sampleCounts = mSampleCounts.getCol(i);
            final double[] sigAllocs = sampleContribs.getCol(i);

            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal == 0)
                continue;

            final double[] fittedCounts = calculateFittedCounts(mSignatures, sigAllocs);
            SigResiduals residuals = calcResiduals(sampleCounts, fittedCounts, sampleTotal);

            double allocTotal = sumVector(sigAllocs);
            double allocPerc = allocTotal / sampleTotal;

            SIG_LOGGER.debug(String.format("sample(%s) alloc(%s perc=%.3f of total=%s) residuals(%.3f total=%s excess=%s)",
                    sampleId, sizeToStr(allocTotal), allocPerc, sizeToStr(sampleTotal),
                    residuals.Percent, sizeToStr(residuals.Total), sizeToStr(residuals.Excess)));

            totalCounts +=sampleTotal;
            totalAlloc += sumVector(sigAllocs);
            totalResiduals += residuals.Total;
            totalResidualsExcess += residuals.Excess;

            List<Integer> sortedSigs = getSortedVectorIndices(sigAllocs, false);
            for (Integer sigIndex : sortedSigs)
            {
                double sigAlloc = sigAllocs[sigIndex];
                double sigPercent = sigAlloc / sampleTotal;
                final String sigName = mSigNames.get(sigIndex);

                writeSampleResults(fitMethod, sampleId, sigName, sigAlloc, sigPercent);

                SIG_LOGGER.trace(String.format("sample(%s) sampleTotal(%.0f) sig(%d) alloc(%.0f perc=%.3f)",
                        sampleId, sampleTotal, sigIndex, sigAlloc, sigPercent));

                if (sigAlloc < 1)
                    break;
            }

            writeSampleResults(fitMethod, sampleId, SIG_UNALLOCATED, residuals.unallocated(), residuals.unallocated()/sampleTotal);
            writeSampleResults(fitMethod, sampleId, SIG_EXCESS, residuals.Excess, residuals.Excess/sampleTotal);

            analyseBucketResiduals(fitMethod, sampleId, sampleCounts, fittedCounts, sigAllocs);
        }

        SIG_LOGGER.info(String.format("method(%s) summary: samples(%s) total(%s) alloc(%.3f %s) residuals(%s excess=%s)",
                fitMethod, sampleCount, sizeToStr(totalCounts), totalAlloc/totalCounts, sizeToStr(totalAlloc),
                sizeToStr(totalResiduals), sizeToStr(totalResidualsExcess)));
    }

    private void analyseBucketResiduals(
            FitMethod fitMethod, final String sampleId, final double[] sampleCounts, final double[] fittedCounts, final double[] sigAllocs)
    {
        for(int bucket = 0; bucket < sampleCounts.length; ++bucket)
        {
            int bucketCount = (int)sampleCounts[bucket];

            if(bucketCount == 0)
                continue;

            double fittedCount = fittedCounts[bucket];
            double fitDiff = bucketCount - fittedCount;
            double fitDiffPerc = abs(fitDiff) / bucketCount;

            if(fitDiffPerc < 0.01)
                continue;

            // bucket has been over-fitted - ie more allocated to signatures than the actual count
            int permittedNoise = calcRangeValue(mNoiseRangeMap, bucketCount, mNoiseProbability, DEFAULT_MIN_NOISE_COUNT, true);

            if(abs(fitDiff) <= permittedNoise)
                continue;

            if(fitDiff < 0)
            {
                // find the sig(s) which contribute the most to this overfit / excess
                for(int s = 0; s < sigAllocs.length; ++s)
                {
                    double sigBucketFit = mSignatures.getCol(s)[bucket] * sigAllocs[s];

                    if(sigBucketFit > abs(fitDiff))
                    {
                        writeBucketFitResult(fitMethod, sampleId, bucket, bucketCount, fittedCount, mSigNames.get(s), sigBucketFit);
                    }
                }
            }
            else
            {
                writeBucketFitResult(fitMethod, sampleId, bucket, bucketCount, fittedCount, "NONE", 0);
            }
        }
    }

    private void intialiseOutputFiles()
    {
        try
        {
            mSampleResultsWriter = createBufferedWriter(formOutputFilename(mOutputDir, mOutputId, "sig_fit_results"), false);

            mSampleResultsWriter.write("SampleId,FitMethod,SigName,Allocation,Percent");
            mSampleResultsWriter.newLine();

            mResidualsWriter = createBufferedWriter(formOutputFilename(mOutputDir, mOutputId, "sig_fit_residuals"), false);
            mResidualsWriter.write("SampleId,FitMethod,Bucket,Count,FitCount,SigName,SigFitCount");
            mResidualsWriter.newLine();

        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to initialise output files: {}", e.toString());
        }
    }

    private void writeSampleResults(
            final FitMethod fitMethod, final String sampleId, final String sigName, double allocTotal, double allocPercent)
    {
        try
        {
            mSampleResultsWriter.write(String.format("%s,%s,%s,%.1f,%.4f", sampleId, fitMethod, sigName, allocTotal, allocPercent));
            mSampleResultsWriter.newLine();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write sample sig results: {}", e.toString());
        }
    }

    private void writeBucketFitResult(
            final FitMethod fitMethod, final String sampleId, int bucket, int bucketCount,
            double fitCount, final String sigName, double sigFitCount)
    {
        try
        {
            mResidualsWriter.write(String.format("%s,%s,%d,%d,%.1f,%s,%.1f",
                    sampleId, fitMethod, bucket, bucketCount, fitCount, sigName, sigFitCount));
            mResidualsWriter.newLine();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write bucket fit results: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        CommonUtils.addCmdLineArgs(options);
        options.addOption(FIT_METHODS, true, "Signatures fit method: NMF, Bucket, LeastSquares");
        BaSampleFitter.addCmdLineArgs(options);

        NmfConfig.addCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        FitAnalyser fitAnalyser = new FitAnalyser(cmd);
        fitAnalyser.run();
    }
}
