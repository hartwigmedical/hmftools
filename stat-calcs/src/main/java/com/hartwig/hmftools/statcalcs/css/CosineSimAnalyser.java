package com.hartwig.hmftools.statcalcs.css;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.sigs.NoiseCalcs.calcPoissonRangeGivenProb;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.statcalcs.common.StatsCommon.LOG_DEBUG;
import static com.hartwig.hmftools.statcalcs.common.StatsCommon.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.statcalcs.common.StatsCommon.STAT_LOGGER;
import static com.hartwig.hmftools.statcalcs.common.StatsCommon.formOutputFilename;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.statcalcs.common.StatsCommon;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CosineSimAnalyser
{
    // config
    private final String mOutputDir;
    private final String mOutputId;

    private final Matrix mSampleCounts;
    private final List<String> mSampleIds;
    private final Map<String,String> mSampleCancerTypes;
    private final Map<String,Integer> mSampleCountsIndex;
    private final double mCssThreshold;
    private final boolean mUseElevated;
    private final Map<Integer,Integer> mRangeMap;

    private final Matrix mReferenceSampleCounts;
    private final List<String> mRefNames;

    private BufferedWriter mWriter;

    private static final String CSS_THRESHOLD = "css_threshold";
    private static final String USE_ELEVATED = "use_elevated";
    private static final String SAMPLE_COUNTS_FILE = "sample_counts_file";
    private static final String REF_COUNTS_FILE = "ref_counts_file";
    private static final String SAMPLE_REF_FILE = "sample_ref_file";

    private static final int MIN_ELEVATED_BUCKETS = 10;
    private static final double ELEVATED_BUCKET_PROBABILITY = 0.1;
    private static final double DEFAULT_CSS_THRESHOLD = 0.8;

    public CosineSimAnalyser(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_FILE_ID);

        mSampleCountsIndex = Maps.newHashMap();

        mCssThreshold = cmd.hasOption(CSS_THRESHOLD) ? Double.parseDouble(cmd.getOptionValue(CSS_THRESHOLD)) : DEFAULT_CSS_THRESHOLD;
        mUseElevated = cmd.hasOption(USE_ELEVATED);
        mRangeMap = Maps.newHashMap();

        mSampleCounts = loadMatrixDataFile(cmd.getOptionValue(SAMPLE_COUNTS_FILE), mSampleCountsIndex, null, false);
        mSampleCounts.cacheTranspose();

        mSampleCancerTypes = Maps.newHashMap();
        mSampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_REF_FILE))
        {
            loadSampleRefDataFile(cmd.getOptionValue(SAMPLE_REF_FILE));
        }
        else
        {
            mSampleIds.addAll(mSampleCountsIndex.keySet());
        }

        mRefNames = Lists.newArrayList();
        if(cmd.hasOption(REF_COUNTS_FILE))
        {
            STAT_LOGGER.info("loading reference data from file({})", cmd.getOptionValue(REF_COUNTS_FILE));

            mReferenceSampleCounts = loadMatrixDataFile(cmd.getOptionValue(REF_COUNTS_FILE), mRefNames, false);
            mReferenceSampleCounts.cacheTranspose();
        }
        else
        {
            mReferenceSampleCounts = null;
        }

        mWriter = null;
    }

    private void loadSampleRefDataFile(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            for(final String line : fileData)
            {
                final String[] items = line.split(",", -1);
                final String sampleId = items[fieldsIndexMap.get("SampleId")];
                final String cancerType = items[fieldsIndexMap.get("CancerType")];
                mSampleIds.add(sampleId);
                mSampleCancerTypes.put(sampleId, cancerType);
            }

            STAT_LOGGER.info("loaded {} samples from file({})", filename);
        }
        catch (IOException e)
        {
            STAT_LOGGER.error("failed to read sample data file({}): {}", filename, e.toString());
        }
    }

    public void run()
    {
        if(mReferenceSampleCounts != null && mReferenceSampleCounts.Rows != mSampleCounts.Rows)
        {
            STAT_LOGGER.error("sampleCounts buckets({}) refCounts buckets({}) mismatch", mReferenceSampleCounts.Rows, mSampleCounts.Rows);
            return;
        }

        STAT_LOGGER.info("running CSS comparison for {} samples", mSampleIds.size());

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            final String sampleId = mSampleIds.get(i);

            if(mSampleCountsIndex.get(sampleId) == null)
                continue;

            final double[] sampleCounts = mSampleCounts.getCol(mSampleCountsIndex.get(sampleId));

            if(mReferenceSampleCounts == null)
            {
                for(int j = i + 1; j < mSampleIds.size(); ++j)
                {
                    final String sampleId2 = mSampleIds.get(j);
                    final double[] sampleCounts2 = mSampleCounts.getCol(mSampleCountsIndex.get(sampleId2));

                    double css = mUseElevated ?
                            calcElevatedCountsCss(sampleCounts, sampleCounts2) : calcCosineSim(sampleCounts, sampleCounts2);

                    if(css >= mCssThreshold)
                    {
                        writeCssResults(sampleId, sampleId2, css);
                    }
                }
            }
            else
            {
                double sampleTotal = sumVector(sampleCounts);
                final String sampleCancerType = mSampleCancerTypes.get(sampleId);

                for(int j = 0; j < mRefNames.size(); ++j)
                {
                    final String refName = mRefNames.get(j);
                    final double[] refCounts = mReferenceSampleCounts.getCol(j);

                    double css;

                    if(sampleCancerType != null && refName.equals(sampleCancerType))
                    {
                        double[] adjustedRefCounts = new double[refCounts.length];
                        copyVector(refCounts, adjustedRefCounts);

                        for(int b = 0; b < refCounts.length; ++b)
                        {
                            adjustedRefCounts[b] = max(adjustedRefCounts[b] - (sampleCounts[b]), 0);
                        }

                        css = calcCosineSim(sampleCounts, adjustedRefCounts);
                    }
                    else
                    {
                        css = calcCosineSim(sampleCounts, refCounts);
                    }

                    if(css >= mCssThreshold)
                    {
                        writeCssResults(sampleId, refName, css);
                    }
                }
            }

            if(i > 0 && (i % 100) == 0)
            {
                STAT_LOGGER.info("processed {} samples", i);
            }
        }

        STAT_LOGGER.info("CSS comparison complete");

        closeBufferedWriter(mWriter);
    }

    private double calcElevatedCountsCss(final double[] sampleCounts1, final double[] sampleCounts2)
    {
        double[] elevCounts1 = extractElevatedCounts(sampleCounts1);
        double[] elevCounts2 = extractElevatedCounts(sampleCounts2);
        int elevatedBuckets = 0;

        for(int b = 0; b < sampleCounts1.length; ++b)
        {
            if(elevCounts1[b] > 0 || elevCounts2[b] > 0)
            {
                elevCounts1[b] = sampleCounts1[b];
                elevCounts2[b] = sampleCounts2[b];
                ++elevatedBuckets;
            }
        }

        if(elevatedBuckets < MIN_ELEVATED_BUCKETS)
            return 0;

        return calcCosineSim(elevCounts1, elevCounts2);
    }

    private double[] extractElevatedCounts(final double[] counts)
    {
        double total = sumVector(counts);
        int bucketCount = counts.length;
        int expectedValue = (int)round(total / bucketCount);

        Integer maxPermittedValue = mRangeMap.get(expectedValue);

        if(maxPermittedValue == null)
        {
            int permittedRange = calcPoissonRangeGivenProb(expectedValue, ELEVATED_BUCKET_PROBABILITY, 0, false);
            maxPermittedValue = expectedValue + permittedRange;
            mRangeMap.put(expectedValue, maxPermittedValue);
        }

        double[] elevatedCounts = new double[bucketCount];

        for(int b = 0; b < bucketCount; ++b)
        {
            if(counts[b] > maxPermittedValue)
                elevatedCounts[b] = counts[b];
        }

        return elevatedCounts;
    }

    private void writeCssResults(
            final String sampleId1, final String sampleId2, double css)
    {
        try
        {
            if(mWriter == null)
            {
                mWriter = createBufferedWriter(formOutputFilename(mOutputDir, mOutputId, "sig_css"), false);

                if(mRefNames.isEmpty())
                    mWriter.write("SampleId1,SampleId2,CSS");
                else
                    mWriter.write("SampleId,RefName,CSS");

                mWriter.newLine();
            }

            mWriter.write(String.format("%s,%s,%.6f", sampleId1, sampleId2, css));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            STAT_LOGGER.error("failed to write CSS results: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        StatsCommon.addCmdLineArgs(options);
        options.addOption(CSS_THRESHOLD, true, "Optional - min CSS to log (default = 0.8)");
        options.addOption(USE_ELEVATED, false, "Optional - only include elevated counts in comparison");
        options.addOption(REF_COUNTS_FILE, true, "Optional - min CSS to log (default = 0.8)");
        options.addOption(SAMPLE_REF_FILE, true, "Optional - sample to ref-type mapping file");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        CosineSimAnalyser cosineSims = new CosineSimAnalyser(cmd);
        cosineSims.run();
    }
}
