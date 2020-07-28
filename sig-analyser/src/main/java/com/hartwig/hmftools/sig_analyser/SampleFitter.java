package com.hartwig.hmftools.sig_analyser;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.GENERIC_INPUT_FILE;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.SigUtils.getNewFile;
import static com.hartwig.hmftools.common.sigs.DataUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.sigs.SigMatrix.writeMatrixData;
import static com.hartwig.hmftools.sig_analyser.nmf.NmfConfig.NMF_EXIT_LEVEL;
import static com.hartwig.hmftools.sig_analyser.nmf.NmfConfig.NMF_MAX_ITERATIONS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.ExpectationMaxFit;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.sig_analyser.buckets.BaSampleFitter;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.sig_analyser.nmf.NmfConfig;
import com.hartwig.hmftools.sig_analyser.nmf.NmfSampleFitter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SampleFitter
{
    private static final String FIT_METHOD = "fit_method";
    private static final String FIT_METHOD_NMF = "NMF";
    private static final String FIT_METHOD_BUCKET = "Bucket";
    private static final String FIT_METHOD_LEAST_SQ = "LeastSquares";
    private static final String FIT_METHOD_EXPECTATION_MAX = "ExpectationMax";
    private static final String SIGNATURES_FILE = "signatures_file";

    private static final Logger LOGGER = LogManager.getLogger(SampleFitter.class);

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(GENERIC_INPUT_FILE, true, "Path to the main input file");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(FIT_METHOD, true, "Signatures fit method: NMF, Bucket, LeastSquares");
        options.addOption(SIGNATURES_FILE, true, "Signature definitions");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        BaSampleFitter.addCmdLineArgs(options);

        NmfConfig.addCmdLineArgs(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        final GenericDataCollection scCollection = GenericDataLoader.loadFile(cmd.getOptionValue(GENERIC_INPUT_FILE));
        final SigMatrix sampleCountsMatrix = DataUtils.createMatrixFromListData(scCollection.getData());

        final GenericDataCollection sigsCollection = GenericDataLoader.loadFile(cmd.getOptionValue(SIGNATURES_FILE));
        final SigMatrix signatures = DataUtils.createMatrixFromListData(sigsCollection.getData());

        int sampleCount = sampleCountsMatrix.Cols;
        int sigCount = signatures.Cols;
        final SigMatrix sampleContribs = new SigMatrix(sigCount, sampleCount);

        final String fitMethod = cmd.getOptionValue(FIT_METHOD);

        LOGGER.info("method({}) fitting {} samples with {} signatures", fitMethod, sampleCount, sigCount);

        if(fitMethod.equals(FIT_METHOD_LEAST_SQ))
        {
            LeastSquaresFit lsqFit = new LeastSquaresFit(signatures.Rows, signatures.Cols);

            for(int i = 0; i < sampleCount; ++i)
            {
                final double[] sampleCounts = sampleCountsMatrix.getCol(i);
                lsqFit.initialise(signatures.getData(), sampleCounts);
                lsqFit.solve();

                final double[] sigAllocs = lsqFit.getContribs();
                sampleContribs.setCol(i, sigAllocs);
            }
        }
        else if(fitMethod.equals(FIT_METHOD_EXPECTATION_MAX))
        {
            for(int i = 0; i < sampleCount; ++i)
            {
                final double[] sampleCounts = sampleCountsMatrix.getCol(i);
                final double[] sigAllocs = ExpectationMaxFit.performFit(sampleCounts, signatures, 0.001, 100);
                sampleContribs.setCol(i, sigAllocs);
            }
        }
        else if(fitMethod.equals(FIT_METHOD_NMF))
        {
            NmfConfig nmfConfig = new NmfConfig(
                    Double.parseDouble(cmd.getOptionValue(NMF_EXIT_LEVEL)),
                    Integer.parseInt(cmd.getOptionValue(NMF_MAX_ITERATIONS, "100")));

            NmfSampleFitter nmfFitter = new NmfSampleFitter(nmfConfig, sampleCountsMatrix, signatures);

            nmfFitter.fitSamples();

            sampleContribs.setData(nmfFitter.getContributions().getData());
        }
        else if(fitMethod.equals(FIT_METHOD_BUCKET))
        {
            BaSampleFitter sampleFitter = new BaSampleFitter(sampleCountsMatrix, signatures, cmd);
            sampleFitter.fitAllSamples();
            sampleContribs.setData(sampleFitter.getContributions().getData());
        }

        if(LOGGER.isDebugEnabled())
        {
            for(int i = 0; i < sampleCount; ++i)
            {
                final double[] sampleCounts = sampleCountsMatrix.getCol(0);
                final double[] sigAllocs = sampleContribs.getCol(i);

                double sampleTotal = sumVector(sampleCounts);
                List<Integer> sortedSigs = getSortedVectorIndices(sigAllocs, false);
                for (Integer sigIndex : sortedSigs)
                {
                    double sigAlloc = sigAllocs[sigIndex];
                    double sigPercent = sigAlloc / sampleTotal;

                    LOGGER.debug(String.format("sample(%s) sampleTotal(%.0f) sig(%d) alloc(%.0f perc=%.3f)",
                            scCollection.getFieldNames().get(i), sampleTotal, sigIndex, sigAlloc, sigPercent));

                    if (sigAlloc < 1)
                        break;
                }
            }
        }

        final String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        final String outputFileId = cmd.getOptionValue(OUTPUT_FILE_ID, "siga");
        final String outputFile = outputFileId + "_sample_contribs.csv";

        try
        {
            BufferedWriter writer = getNewFile(outputDir, outputFile);
            writeMatrixData(writer, scCollection.getFieldNames(), sampleContribs, false);
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write sample contributions: {}", e.toString());
            return;
        }

        LOGGER.info("sample signature contributions written");

    }
}
