package com.hartwig.hmftools.sig_analyser.fitter;

import static com.hartwig.hmftools.common.sigs.DataUtils.round;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_UNALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleListFile;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.loadSampleMatrixCounts;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sig_analyser.loaders.SigSnvLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SampleFitter
{
    // config
    private static final String SIGNATURES_FILE = "signatures_file";
    private static final String SAMPLE_IDS = "sample";
    
    private final String mSampleIdsConfig;
    private final List<String> mSampleIdList;
    private final String mSnvCountsFile;
    private final String mSignaturesFile;
    private final String mOutputDir;

    private SigMatrix mSampleCountsMatrix;
    private final List<String> mSignatureNames;
    private SigMatrix mSignatures;
    private DatabaseAccess mDbAccess;
    private BufferedWriter mFitWriter;

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
        mFitWriter = null;
        mDbAccess = createDatabaseAccess(cmd);
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

        final GenericDataCollection sigsCollection = GenericDataLoader.loadFile(mSignaturesFile);
        mSignatures = DataUtils.createMatrixFromListData(sigsCollection.getData());
        mSignatureNames.addAll(sigsCollection.getFieldNames());

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
        snvLoader.loadData(mDbAccess);

        final String filename = mSampleIdList.size() == 1 ? mSampleIdList.get(0) + ".sig.snv_counts.csv" : "SIG_SNV_COUNTS.csv";
        snvLoader.writeSampleCounts(mOutputDir + filename);

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

        final String filename = mOutputDir + "SIG_SNV_FIT.csv";

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
                SignatureAllocationFile.write(SignatureAllocationFile.generateFilename(mOutputDir, sampleId), sigAllocations);
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
        options.addOption(SAMPLE_IDS, true, "Either a single sampleId or a file with a list");
        options.addOption(SAMPLE_COUNTS_FILE, true, "Path to the main input file");
        options.addOption(SIGNATURES_FILE, true, "Signature definitions");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        addDatabaseCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SampleFitter sampleFitter = new SampleFitter(cmd);
        sampleFitter.run();
    }
}
