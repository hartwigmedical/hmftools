package com.hartwig.hmftools.sigs.fitter;

import static com.hartwig.hmftools.common.sigs.DataUtils.round;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_MISALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.utils.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sigs.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;
import static com.hartwig.hmftools.sigs.common.CommonUtils.loadSampleMatrixCounts;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.sigs.loaders.SigSnvLoader;

import org.jetbrains.annotations.NotNull;

public class SignatureFitter
{
    private final SigsConfig mConfig;
    private final VariantFilters mFilters;
    private Matrix mSampleCountsMatrix;
    private final List<String> mSignatureNames;
    private Matrix mSignatures;

    private final boolean mFitToTotal;

    private final SigSnvLoader mSnvLoader;
    private BufferedWriter mFitWriter;

    private static final int ALLOC_ROUND = 3;
    private static final int PERC_ROUND = 5;

    public SignatureFitter(final ConfigBuilder configBuilder)
    {
        mConfig = new SigsConfig(configBuilder);

        mFilters = new VariantFilters(configBuilder);
        
        mSignatureNames = Lists.newArrayList();

        mSampleCountsMatrix = null;
        mSignatures = null;

        mFitWriter = null;
        mFitToTotal = true;

        mSnvLoader = new SigSnvLoader(mFilters);

        if(mConfig.PositionBucketSize > 0)
            mSnvLoader.initialisePositionFrequencies(null, Lists.newArrayList(mConfig.PositionBucketSize));
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
        if(mConfig.SampleIdList.isEmpty())
        {
            SIG_LOGGER.error("no sample ID(s) configured");
            return false;
        }

        if(mConfig.VcfFile == null && mConfig.SnvCountsFile == null)
        {
            SIG_LOGGER.error("no sample VCF or sample count file configured");
            return false;
        }

        if(mConfig.SignaturesFile == null)
        {
            SIG_LOGGER.info("using default SNV Cosmic signatures");
            mSignatures = SnvSigUtils.loadSnvSignatures(mSignatureNames);
        }
        else
        {
            mSignatures = loadMatrixDataFile(mConfig.SignaturesFile, mSignatureNames, false);
        }

        if(mConfig.SnvCountsFile != null)
        {
            mSampleCountsMatrix = loadSampleMatrixCounts(mConfig.SnvCountsFile, mConfig.SampleIdList);
        }

        initialiseOutputFiles();
        return true;
    }

    private void performFit()
    {
        int sampleCount = mConfig.SampleIdList.size();
        int sigCount = mSignatures.Cols;

        Matrix sampleContribs = new Matrix(sigCount, sampleCount);

        SIG_LOGGER.info("fitting sample({}) with {} signatures", sampleCount, sigCount);

        LeastSquaresFit lsqFit = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);

        for(int i = 0; i < sampleCount; ++i)
        {
            String sampleId = mConfig.SampleIdList.get(i);
            double[] sampleCounts = getSampleCounts(sampleId, i);

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

        mSnvLoader.loadData(mConfig.VcfFile, false);

        if(mConfig.SampleIdList.size() == 1)
        {
            final String filename = formOutputFilename(mConfig.OutputDir, mConfig.OutputId, sampleId + ".sig.snv_counts");

            mSnvLoader.writeSampleCounts(filename);

            if(mConfig.PositionBucketSize > 0)
            {
                PositionFreqBuilder posFreqBuilder = new PositionFreqBuilder(
                        mConfig.RefGenVersion, mConfig.OutputDir, mConfig.OutputId, mConfig.PositionBucketSize);

                final PositionFrequencies samplePosFrequencies = mSnvLoader.getPositionFrequencies().get(0);
                posFreqBuilder.writeSampleCounts(sampleId, samplePosFrequencies.getChrPosBucketFrequencies(), mConfig.WritePosFreqCoords);
            }
        }

        return mSnvLoader.getSampleBucketCounts().getCol(0);
    }

    private void processSampleResults(final String sampleId, final double[] sampleCounts, double sampleTotal, final double[] sigAllocs)
    {
        List<SignatureAllocation> sigAllocations = Lists.newArrayList();

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
            if(sigAlloc < mConfig.MinFit || allocPerc < mConfig.MinFitPerc)
                continue;

            String sigName = mSignatureNames.get(s);

            sigAllocations.add(ImmutableSignatureAllocation.builder()
                    .signature(sigName)
                    .allocation(round(sigAlloc, ALLOC_ROUND))
                    .percent(round(allocPerc, PERC_ROUND))
                    .build());
        }

        // also report residuals
        double[] fittedCounts = calculateFittedCounts(mSignatures, sigAllocs);
        SigResiduals residuals = calcResiduals(sampleCounts, fittedCounts, sampleTotal);
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
        if(mConfig.SampleIdList.size() <= 1)
            return;

        final String filename = formOutputFilename(mConfig.OutputDir, mConfig.OutputId, "sig_snv_fit");

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
        boolean isSingleSample = mConfig.SampleIdList.size() == 1;

        try
        {
            if(isSingleSample)
            {
                String fileId = mConfig.OutputId != null ? sampleId + "." + mConfig.OutputId : sampleId;
                SignatureAllocationFile.write(SignatureAllocationFile.generateFilename(mConfig.OutputDir, fileId), sigAllocations);
            }
            else
            {
                for(SignatureAllocation sigAlloc : sigAllocations)
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        SigsConfig.registerConfig(configBuilder);
        VariantFilters.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SignatureFitter signatureFitter = new SignatureFitter(configBuilder);
        signatureFitter.run();
    }
}
