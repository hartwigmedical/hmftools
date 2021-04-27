package com.hartwig.hmftools.sigs.fitter;

import static com.hartwig.hmftools.common.sigs.DataUtils.round;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.PositionFrequencies.DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_UNALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.VectorUtils.vectorMultiply;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.sigs.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SAMPLE_IDS;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIGNATURES_FILE;
import static com.hartwig.hmftools.sigs.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sigs.common.CommonUtils.formOutputFilename;
import static com.hartwig.hmftools.sigs.common.CommonUtils.loadSampleListFile;
import static com.hartwig.hmftools.sigs.common.CommonUtils.loadSampleMatrixCounts;
import static com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder.MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder.POSITION_BUCKET_SIZE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sigs.common.CommonUtils;
import com.hartwig.hmftools.sigs.loaders.PositionFreqBuilder;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.sigs.loaders.SigSnvLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SampleFitter
{
    private final String mSampleIdsConfig;
    private final List<String> mSampleIdList;
    private final String mSnvCountsFile;
    private final String mSignaturesFile;
    private final String mOutputDir;
    private final String mOutputId;

    private Matrix mSampleCountsMatrix;
    private final List<String> mSignatureNames;
    private Matrix mSignatures;

    private final int mPositionBucketSize;
    private final int mMaxSampleCount;
    private final boolean mFitToTotal;
    private final boolean mWritePosFreqCoords;

    private final DatabaseAccess mDbAccess;
    private final String mVcfFile;
    private final SigSnvLoader mSnvLoader;
    private boolean mUploadToDb;
    private BufferedWriter mFitWriter;

    private static double MIN_ALLOCATION = 1;
    private static double MIN_ALLOCATION_PERC = 0.005;
    private static final int ALLOC_ROUND = 3;
    private static final int PERC_ROUND = 5;

    private static final String UPLOAD_TO_DB = "upload_to_db";
    private static final String SOMATIC_VCF_FILE = "somatic_vcf_file";
    private static final String MIN_ALLOC = "min_alloc";
    private static final String MIN_ALLOC_PERC = "min_alloc_perc";
    private static final String FIT_TO_TOTAL = "fit_to_total";
    private static final String WRITE_POS_COORDS = "write_pos_freq_coords";

    public SampleFitter(final CommandLine cmd)
    {
        mSnvCountsFile = cmd.getOptionValue(SAMPLE_COUNTS_FILE);
        mSignaturesFile = cmd.getOptionValue(SIGNATURES_FILE);
        mSampleIdsConfig = cmd.getOptionValue(SAMPLE_IDS);
        mSampleIdList = Lists.newArrayList();
        mSignatureNames = Lists.newArrayList();

        mSampleCountsMatrix = null;
        mSignatures = null;

        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mFitWriter = null;
        mDbAccess = createDatabaseAccess(cmd);
        mVcfFile = cmd.getOptionValue(SOMATIC_VCF_FILE);
        mUploadToDb = cmd.hasOption(UPLOAD_TO_DB);
        mFitToTotal = Boolean.parseBoolean(cmd.getOptionValue(FIT_TO_TOTAL, "true"));
        mWritePosFreqCoords = cmd.hasOption(WRITE_POS_COORDS);

        mPositionBucketSize = Integer.parseInt(cmd.getOptionValue(POSITION_BUCKET_SIZE, "0"));
        mMaxSampleCount = Integer.parseInt(cmd.getOptionValue(MAX_SAMPLE_COUNT, String.valueOf(DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT)));

        MIN_ALLOCATION = Double.parseDouble(cmd.getOptionValue(MIN_ALLOC, "0.01"));
        MIN_ALLOCATION_PERC = Double.parseDouble(cmd.getOptionValue(MIN_ALLOC_PERC, "0.0005"));

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

        if(mSampleIdsConfig == null && mSnvCountsFile == null)
        {
            SIG_LOGGER.error("missing sampleIds config and sample count file");
            return false;
        }

        mSignatures = loadMatrixDataFile(mSignaturesFile, mSignatureNames);

        if(mVcfFile != null)
        {
            mSampleIdList.add(mSampleIdsConfig);
        }
        else
        {
            if(mSnvCountsFile != null)
            {
                mSampleCountsMatrix = loadSampleMatrixCounts(mSnvCountsFile, mSampleIdList);
            }
            else
            {
                if(mDbAccess == null)
                {
                    SIG_LOGGER.error("missing DB connection when no sample counts file configured");
                    return false;
                }

                // load from file or delimitered list
                if(mSampleIdsConfig.contains(".csv"))
                {
                    // load from file
                    mSampleIdList.addAll(loadSampleListFile(mSampleIdsConfig));
                }
                else
                {
                    mSampleIdList.add(mSampleIdsConfig);
                }
            }
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

            if(sampleTotal == 0)
                continue;

            lsqFit.initialise(mSignatures.getData(), sampleCounts);
            lsqFit.solve();

            final double[] sigAllocs = lsqFit.getContribs();
            sampleContribs.setCol(i, sigAllocs);

            processSampleResults(sampleId, sampleCounts, sampleTotal, sigAllocs);

            if(i > 0 && (i % 100) == 0)
            {
                SIG_LOGGER.info("processed {} samples", i);
            }
        }
    }

    private final double[] getSampleCounts(final String sampleId, int sampleIndex)
    {
        if(mSampleCountsMatrix != null)
            return mSampleCountsMatrix.getCol(sampleIndex);

        mSnvLoader.setSampleIds(Lists.newArrayList(sampleId));

        mSnvLoader.loadData(mDbAccess, mVcfFile, false);

        if(mSampleIdList.size() == 1)
        {
            final String filename = formOutputFilename(mOutputDir, mOutputId, sampleId + ".sig.snv_counts");

            mSnvLoader.writeSampleCounts(filename);

            if(mPositionBucketSize > 0)
            {
                PositionFreqBuilder posFreqBuilder =
                        new PositionFreqBuilder(mOutputDir, mOutputId, mPositionBucketSize);

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
            if(sigAlloc < MIN_ALLOCATION || allocPerc < MIN_ALLOCATION_PERC)
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
                .signature(SIG_UNALLOCATED)
                .allocation(round(residualsUnalloc,ALLOC_ROUND))
                .percent(round(residualsUnalloc/sampleTotal, PERC_ROUND))
                .build());

        sigAllocations.add(ImmutableSignatureAllocation.builder()
                .signature(SIG_EXCESS)
                .allocation(round(residuals.Excess,ALLOC_ROUND))
                .percent(round(residuals.Excess/sampleTotal, PERC_ROUND))
                .build());

        writeSigAllocations(sampleId, sigAllocations);

        if(mDbAccess != null && mUploadToDb)
        {
            SIG_LOGGER.debug("sample({}) writing {} allocations to database", sampleId, sigAllocations.size());
            mDbAccess.writeSignatures(sampleId, sigAllocations);
        }
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
            SIG_LOGGER.error("error writing to outputFile({}): {}", filename, e.toString());
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
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        CommonUtils.addCmdLineArgs(options);
        addDatabaseCmdLineArgs(options);
        options.addOption(POSITION_BUCKET_SIZE, true, "Position bucket size");
        options.addOption(MAX_SAMPLE_COUNT, true, "Max sample SNV count for position frequencies, default = 20K");
        options.addOption(UPLOAD_TO_DB, false, "Upload results to database if present");
        options.addOption(SOMATIC_VCF_FILE, true, "Somatic variant VCF file");
        options.addOption(MIN_ALLOC_PERC, true, "Min signature allocation as percentage (default=0.5%)");
        options.addOption(MIN_ALLOC, true, "Min signature allocation (default=1)");
        options.addOption(WRITE_POS_COORDS, false, "Include coordinates with pos-frequency counts");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SampleFitter sampleFitter = new SampleFitter(cmd);
        sampleFitter.run();
    }
}
