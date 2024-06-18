package com.hartwig.hmftools.sigs.fitter;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.sigs.DataUtils.round;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_MISALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sigs.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIGNATURES_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;
import static com.hartwig.hmftools.sigs.common.CommonUtils.loadSampleMatrixCounts;
import static com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder.MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder.POSITION_BUCKET_SIZE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.sigs.loaders.SigSnvLoader;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SampleFitter
{
    private final RefGenomeVersion mRefGenomeVersion;
    private final List<String> mSampleIdList;
    private final String mSnvCountsFile;
    private final String mSignaturesFile;
    private final String mOutputDir;
    private final String mOutputId;

    private Matrix mSampleCountsMatrix;
    private final List<String> mSignatureNames;
    private Matrix mSignatures;

    private final int mPositionBucketSize;
    private final boolean mFitToTotal;
    private final boolean mWritePosFreqCoords;

    private final String mVcfFile;
    private final SigSnvLoader mSnvLoader;
    private BufferedWriter mFitWriter;

    private final double mMinFit;
    private final double mMinFitPerc;

    private static final double DEFAULT_MIN_ALLOCATION = 0.01;
    private static final double DEFAULT_MIN_ALLOCATION_PERC = 0.0005;
    private static final int ALLOC_ROUND = 3;
    private static final int PERC_ROUND = 5;

    private static final String SOMATIC_VCF_FILE = "somatic_vcf_file";
    private static final String MIN_ALLOC = "min_alloc";
    private static final String MIN_ALLOC_PERC = "min_alloc_perc";
    private static final String FIT_TO_TOTAL = "fit_to_total";
    private static final String WRITE_POS_COORDS = "write_pos_freq_coords";

    public SampleFitter(final ConfigBuilder configBuilder)
    {
        mSnvCountsFile = configBuilder.getValue(SAMPLE_COUNTS_FILE);
        mSignaturesFile = configBuilder.getValue(SIGNATURES_FILE);
        mVcfFile = configBuilder.getValue(SOMATIC_VCF_FILE);

        mSampleIdList = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            mSampleIdList.add(configBuilder.getValue(SAMPLE));
        }
        else if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            mSampleIdList.addAll(loadSampleIdsFile(configBuilder));
        }

        mSignatureNames = Lists.newArrayList();

        mSampleCountsMatrix = null;
        mSignatures = null;

        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);

        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mFitWriter = null;
        mFitToTotal = true; // configBuilder.hasFlag(FIT_TO_TOTAL);
        mWritePosFreqCoords = configBuilder.hasFlag(WRITE_POS_COORDS);

        mPositionBucketSize = configBuilder.getInteger(POSITION_BUCKET_SIZE);
        mMinFit = configBuilder.getDecimal(MIN_ALLOC);
        mMinFitPerc = configBuilder.getDecimal(MIN_ALLOC_PERC);

        mSnvLoader = new SigSnvLoader(null);

        if(mPositionBucketSize > 0)
            mSnvLoader.initialisePositionFrequencies(null, Lists.newArrayList(mPositionBucketSize));
    }

    public void run()
    {
        SIG_LOGGER.info("running sample signature fit");

        if(!initialise())
            return;

        performFit();

        SIG_LOGGER.info("sample signature fit complete");

        closeBufferedWriter(mFitWriter);
    }
    
    private boolean initialise()
    {
        if(mSignaturesFile == null)
        {
            SIG_LOGGER.error("missing mSignatures file");
            return false;
        }

        if(mSampleIdList.isEmpty())
        {
            SIG_LOGGER.error("no sample ID(s) configured");
            return false;
        }

        if(mVcfFile == null && mSnvCountsFile == null)
        {
            SIG_LOGGER.error("no sample VCF or sample count file configured");
            return false;
        }

        mSignatures = loadMatrixDataFile(mSignaturesFile, mSignatureNames, false);

        if(mSnvCountsFile != null)
        {
            mSampleCountsMatrix = loadSampleMatrixCounts(mSnvCountsFile, mSampleIdList);
        }

        initialiseOutputFiles();
        return true;
    }

    private void performFit()
    {
        int sampleCount = mSampleIdList.size();
        int sigCount = mSignatures.Cols;

        final Matrix sampleContribs = new Matrix(sigCount, sampleCount);

        SIG_LOGGER.info("fitting sample({}) with {} signatures", sampleCount, sigCount);

        LeastSquaresFit lsqFit = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);

        for(int i = 0; i < sampleCount; ++i)
        {
            final String sampleId = mSampleIdList.get(i);
            final double[] sampleCounts = getSampleCounts(sampleId, i);

            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal > 0)
            {
                lsqFit.initialise(mSignatures.getData(), sampleCounts);
                lsqFit.solve();

                final double[] sigAllocs = lsqFit.getContribs();
                sampleContribs.setCol(i, sigAllocs);

                processSampleResults(sampleId, sampleCounts, sampleTotal, sigAllocs);
            }
            else
            {
                List<SignatureAllocation> emptySigAllocations = Lists.newArrayList();
                writeSampleSigResults(sampleId, emptySigAllocations);
            }

            if(i > 0 && (i % 100) == 0)
            {
                SIG_LOGGER.info("processed {} samples", i);
            }
        }
    }

    private double[] getSampleCounts(final String sampleId, int sampleIndex)
    {
        if(mSampleCountsMatrix != null)
            return mSampleCountsMatrix.getCol(sampleIndex);

        mSnvLoader.setSampleIds(Lists.newArrayList(sampleId));

        mSnvLoader.loadData(null, mVcfFile, false);

        if(mSampleIdList.size() == 1)
        {
            final String filename = formOutputFilename(mOutputDir, mOutputId, sampleId + ".sig.snv_counts");

            mSnvLoader.writeSampleCounts(filename);

            if(mPositionBucketSize > 0)
            {
                PositionFreqBuilder posFreqBuilder = new PositionFreqBuilder(mRefGenomeVersion, mOutputDir, mOutputId, mPositionBucketSize);

                final PositionFrequencies samplePosFrequencies = mSnvLoader.getPositionFrequencies().get(0);
                posFreqBuilder.writeSampleCounts(sampleId, samplePosFrequencies.getChrPosBucketFrequencies(), mWritePosFreqCoords);
            }
        }

        return mSnvLoader.getSampleBucketCounts().getCol(0);
    }

    private void processSampleResults(final String sampleId, final double[] sampleCounts, double sampleTotal, final double[] sigAllocs)
    {
        final List<SignatureAllocation> sigAllocations = Lists.newArrayList();

        double allocTotal = 0;

        if(mFitToTotal)
        {
            double fitTotal = sumVector(sigAllocs);
            double adjustFactor = sampleTotal / fitTotal;
            vectorMultiply(sigAllocs, adjustFactor);
        }

        for(int s = 0; s < sigAllocs.length; ++s)
        {
            double sigAlloc = sigAllocs[s];
            double allocPerc = sigAlloc / sampleTotal;
            allocTotal += sigAlloc;

            // avoid storing tiny or zero sig allocations
            if(sigAlloc < mMinFit || allocPerc < mMinFitPerc)
                continue;

            final String sigName = mSignatureNames.get(s);

            sigAllocations.add(ImmutableSignatureAllocation.builder()
                    .signature(sigName)
                    .allocation(round(sigAlloc, ALLOC_ROUND))
                    .percent(round(allocPerc, PERC_ROUND))
                    .build());
        }

        // also report residuals
        final double[] fittedCounts = calculateFittedCounts(mSignatures, sigAllocs);
        final SigResiduals residuals = calcResiduals(sampleCounts, fittedCounts, sampleTotal);
        double residualsUnalloc = residuals.Total - residuals.Excess;

        SIG_LOGGER.debug(String.format("sample(%s) alloc(%s perc=%.3f of total=%s) residuals(%.3f total=%s excess=%s)",
                sampleId, sizeToStr(allocTotal), allocTotal/sampleTotal, sizeToStr(sampleTotal),
                residuals.Percent, sizeToStr(residuals.Total), sizeToStr(residuals.Excess)));

        sigAllocations.add(ImmutableSignatureAllocation.builder()
                .signature(SIG_MISALLOCATED)
                .allocation(round(residualsUnalloc,ALLOC_ROUND))
                .percent(round(residualsUnalloc/sampleTotal, PERC_ROUND))
                .build());

        writeSampleSigResults(sampleId, sigAllocations);
    }

    private void writeSampleSigResults(final String sampleId, final List<SignatureAllocation> sigAllocations)
    {
        writeSigAllocations(sampleId, sigAllocations);
    }

    private void initialiseOutputFiles()
    {
        if(mSampleIdList.size() <= 1)
            return;

        final String filename = formOutputFilename(mOutputDir, mOutputId, "sig_snv_fit");

        try
        {
            mFitWriter = createBufferedWriter(filename, false);
            mFitWriter.write("SampleId,Signature,Allocation,Percent");
            mFitWriter.newLine();
        }
        catch (final IOException e)
        {
            SIG_LOGGER.error("error writing to file({}): {}", filename, e.toString());
            System.exit(1);
        }
    }

    private void writeSigAllocations(final String sampleId, final List<SignatureAllocation> sigAllocations)
    {
        boolean isSingleSample = mSampleIdList.size() == 1;

        try
        {
            if(isSingleSample)
            {
                final String fileId = mOutputId != null ? sampleId + "." + mOutputId : sampleId;
                SignatureAllocationFile.write(SignatureAllocationFile.generateFilename(mOutputDir, fileId), sigAllocations);
            }
            else
            {
                for(final SignatureAllocation sigAlloc : sigAllocations)
                {
                    mFitWriter.write(String.format("%s,%s,%g,%g",
                            sampleId, sigAlloc.signature(), sigAlloc.allocation(), sigAlloc.percent()));
                    mFitWriter.newLine();
                }
            }
        }
        catch (final IOException e)
        {
            SIG_LOGGER.error("error writing to fit output: {}", e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        // either a VCF for a single sample or a sample matrix file
        configBuilder.addPath(SOMATIC_VCF_FILE, false, "Somatic variant VCF file");
        configBuilder.addPath(SAMPLE_COUNTS_FILE, false, "Path to the main input file");

        configBuilder.addPath(SIGNATURES_FILE, false, "Signature definitions");
        addRefGenomeVersion(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.addInteger(POSITION_BUCKET_SIZE, "Position bucket size", 0);

        configBuilder.addDecimal(MIN_ALLOC_PERC, "Min signature allocation as percentage", DEFAULT_MIN_ALLOCATION_PERC);
        configBuilder.addDecimal(MIN_ALLOC, "Min signature allocation", DEFAULT_MIN_ALLOCATION);
        configBuilder.addFlag(WRITE_POS_COORDS, "Include coordinates with pos-frequency counts");

        configBuilder.checkAndParseCommandLine(args);

        SampleFitter sampleFitter = new SampleFitter(configBuilder);
        sampleFitter.run();
    }
}
