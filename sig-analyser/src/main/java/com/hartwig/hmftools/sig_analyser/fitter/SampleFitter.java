package com.hartwig.hmftools.sig_analyser.fitter;

import static com.hartwig.hmftools.common.sigs.DataUtils.round;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_UNALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SAMPLE_IDS;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIGNATURES_FILE;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.formOutputFilename;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleListFile;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleMatrixCounts;
import static com.hartwig.hmftools.sig_analyser.loaders.PositionFreqBuilder.DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.sig_analyser.loaders.PositionFreqBuilder.MAX_SAMPLE_COUNT;
import static com.hartwig.hmftools.sig_analyser.loaders.PositionFreqBuilder.POSITION_BUCKET_SIZE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sig_analyser.common.CommonUtils;
import com.hartwig.hmftools.sig_analyser.loaders.PositionFreqBuilder;
import com.hartwig.hmftools.sig_analyser.loaders.PositionFrequencies;
import com.hartwig.hmftools.sig_analyser.loaders.SigSnvLoader;

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

    private SigMatrix mSampleCountsMatrix;
    private final List<String> mSignatureNames;
    private SigMatrix mSignatures;

    private final int mPositionBucketSize;
    private final int mMaxSampleCount;

    private DatabaseAccess mDbAccess;
    private BufferedWriter mFitWriter;

    private static final double MIN_ALLOCATION = 1;
    private static final double MIN_ALLOCATION_PERC = 0.005;

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

        mPositionBucketSize = Integer.parseInt(cmd.getOptionValue(POSITION_BUCKET_SIZE, "0"));
        mMaxSampleCount = Integer.parseInt(cmd.getOptionValue(MAX_SAMPLE_COUNT, String.valueOf(DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT)));
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

        if(mSnvCountsFile != null)
        {
            mSampleCountsMatrix = loadSampleMatrixCounts(mSnvCountsFile, mSampleIdList);
        }
        else
        {
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

            loadSnvCountsData();
        }

        initialiseOutputFiles();
        return true;
    }

    private void loadSnvCountsData()
    {
        if(mDbAccess == null)
            return;

        SigSnvLoader snvLoader = new SigSnvLoader(null, mSampleIdList);

        if(mPositionBucketSize > 0)
        {
            snvLoader.initialisePositionFrequencies(null, Lists.newArrayList(mPositionBucketSize));
        }

        snvLoader.loadData(mDbAccess, false);

        final String filename = mSampleIdList.size() == 1 ?
                formOutputFilename(mOutputDir, mOutputId, mSampleIdList.get(0) + ".sig.snv_counts")
                : formOutputFilename(mOutputDir, mOutputId, "sig_snv_counts");

        snvLoader.writeSampleCounts(filename);

        if(mPositionBucketSize > 0 && mSampleIdList.size() == 1)
        {
            PositionFreqBuilder posFreqBuilder =
                    new PositionFreqBuilder(mOutputDir, mPositionBucketSize, mMaxSampleCount);

            final PositionFrequencies samplePosFrequencies = snvLoader.getPositionFrequencies().get(0);
            posFreqBuilder.writeSampleCounts(mSampleIdList.get(0), samplePosFrequencies.getChrPosBucketFrequencies());
        }

        mSampleCountsMatrix = snvLoader.getSampleBucketCounts();
    }

    private void performFit()
    {
        int sampleCount = mSampleCountsMatrix.Cols;
        int sigCount = mSignatures.Cols;
        
        final SigMatrix sampleContribs = new SigMatrix(sigCount, sampleCount);

        SIG_LOGGER.info("fitting sample({}) with {} signatures", sampleCount, sigCount);

        LeastSquaresFit lsqFit = new LeastSquaresFit(mSignatures.Rows, mSignatures.Cols);

        for(int i = 0; i < sampleCount; ++i)
        {
            final double[] sampleCounts = mSampleCountsMatrix.getCol(i);

            double sampleTotal = sumVector(sampleCounts);

            if(sampleTotal == 0)
                continue;

            lsqFit.initialise(mSignatures.getData(), sampleCounts);
            lsqFit.solve();

            final double[] sigAllocs = lsqFit.getContribs();
            sampleContribs.setCol(i, sigAllocs);

            processSampleResults(mSampleIdList.get(i), sampleCounts, sampleTotal, sigAllocs);
        }
    }

    private void processSampleResults(final String sampleId, final double[] sampleCounts, double sampleTotal, final double[] sigAllocs)
    {
        final List<SignatureAllocation> sigAllocations = Lists.newArrayList();

        double allocTotal = 0;

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
                    .allocation(round(sigAlloc,1))
                    .percent(round(allocPerc, 3))
                    .build());
        }

        // also report residuals
        final double[] fittedCounts = calculateFittedCounts(mSignatures, sigAllocs);
        final SigResiduals residuals = calcResiduals(sampleCounts, fittedCounts, sampleTotal);
        double residualsUnalloc = residuals.Total - residuals.Excess;

        sigAllocations.add(ImmutableSignatureAllocation.builder()
                .signature(SIG_UNALLOCATED)
                .allocation(round(residualsUnalloc,1))
                .percent(round(residualsUnalloc/sampleTotal, 3))
                .build());

        sigAllocations.add(ImmutableSignatureAllocation.builder()
                .signature(SIG_EXCESS)
                .allocation(round(residuals.Excess,1))
                .percent(round(residuals.Excess/sampleTotal, 3))
                .build());

        writeSigAllocations(sampleId, sigAllocations);

        if(mDbAccess != null)
            mDbAccess.writeSignatures(sampleId, sigAllocations);

        SIG_LOGGER.debug(String.format("sample(%s) alloc(%s perc=%.3f of total=%s) residuals(%.3f total=%s excess=%s)",
                sampleId, sizeToStr(allocTotal), allocTotal/sampleTotal, sizeToStr(sampleTotal),
                residuals.Percent, sizeToStr(residuals.Total), sizeToStr(residuals.Excess)));
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
                    mFitWriter.write(String.format("%s,%s,%.1f,%.3f",
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
        options.addOption(MAX_SAMPLE_COUNT, true, "Max sample SNV count, default = 20K");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SampleFitter sampleFitter = new SampleFitter(cmd);
        sampleFitter.run();
    }
}
