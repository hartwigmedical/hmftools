package com.hartwig.hmftools.sig_analyser.analysers;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sigs.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.DataUtils.sizeToStr;
import static com.hartwig.hmftools.common.sigs.SigMatrix.writeMatrixData;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_EXCESS;
import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_UNALLOCATED;
import static com.hartwig.hmftools.common.sigs.SigUtils.calcResiduals;
import static com.hartwig.hmftools.common.sigs.SigUtils.calculateFittedCounts;
import static com.hartwig.hmftools.common.sigs.VectorUtils.getSortedVectorIndices;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.LOG_DEBUG;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SAMPLE_COUNTS_FILE;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIG_LOGGER;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.calcRangeValue;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.ExpectationMaxFit;
import com.hartwig.hmftools.common.sigs.LeastSquaresFit;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.sigs.SigResiduals;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.sig_analyser.buckets.BaSampleFitter;
import com.hartwig.hmftools.sig_analyser.fitter.FitMethod;
import com.hartwig.hmftools.sig_analyser.nmf.NmfConfig;
import com.hartwig.hmftools.sig_analyser.nmf.NmfSampleFitter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CosineSimilarities
{
    // config
    private final String mOutputDir;

    private final SigMatrix mSampleCounts;
    private final List<String> mSampleIds;
    private final Map<String,Integer> mSampleCountsIndex;
    private final double mCssThreshold;

    private BufferedWriter mWriter;

    private static final String SAMPLE_IDS = "samples";
    private static final String CSS_THRESHOLD = "css_threshold";

    public CosineSimilarities(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);

        mSampleCountsIndex = Maps.newHashMap();

        mCssThreshold = Double.parseDouble(cmd.getOptionValue(CSS_THRESHOLD, "0.8"));

        final GenericDataCollection collection = GenericDataLoader.loadFile(cmd.getOptionValue(SAMPLE_COUNTS_FILE));
        mSampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        mSampleCounts.cacheTranspose();

        for(int s = 0; s < collection.getFieldNames().size(); ++s)
        {
            final String sampleId = collection.getFieldNames().get(s);
            mSampleCountsIndex.put(sampleId, s);
        }

        if(cmd.hasOption(SAMPLE_IDS))
        {
            mSampleIds = Arrays.stream(cmd.getOptionValue(SAMPLE_IDS).split(";", -1)).collect(Collectors.toList());
        }
        else
        {
            mSampleIds = collection.getFieldNames();
        }

        mWriter = null;
    }

    public void run()
    {
        SIG_LOGGER.info("running CSS comparison for {} samples", mSampleIds.size());

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            final String sampleId1 = mSampleIds.get(i);
            final double[] sampleCounts1 = mSampleCounts.getCol(mSampleCountsIndex.get(sampleId1));

            for(int j = i + 1; j < mSampleIds.size(); ++j)
            {
                final String sampleId2 = mSampleIds.get(j);
                final double[] sampleCounts2 = mSampleCounts.getCol(mSampleCountsIndex.get(sampleId2));

                double css = calcCosineSim(sampleCounts1, sampleCounts2);

                if(css >= mCssThreshold)
                {
                    writeCssResults(sampleId1, sampleId2, css);
                }
            }

            if(i > 0 && (i % 100) == 0)
            {
                SIG_LOGGER.info("processed {} samples", i);
            }
        }

        SIG_LOGGER.info("CSS comparison complete");

        closeBufferedWriter(mWriter);
    }

    private void writeCssResults(
            final String sampleId1, final String sampleId2, double css)
    {
        try
        {
            if(mWriter == null)
            {
                mWriter = createBufferedWriter(mOutputDir + "SIG_CSS_RESULTS.csv", false);

                mWriter.write("SampleId1,SampleId2,CSS");
                mWriter.newLine();
            }

            mWriter.write(String.format("%s,%s,%.6f", sampleId1, sampleId2, css));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SIG_LOGGER.error("failed to write CSS results: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        options.addOption(SAMPLE_COUNTS_FILE, true, "Sample bucket counts");
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(SAMPLE_IDS, true, "Optional - list of sampleIds, separated by ';");
        options.addOption(CSS_THRESHOLD, true, "Optional - min CSS to log (default = 0.8)");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        CosineSimilarities cosineSims = new CosineSimilarities(cmd);
        cosineSims.run();
    }
}
