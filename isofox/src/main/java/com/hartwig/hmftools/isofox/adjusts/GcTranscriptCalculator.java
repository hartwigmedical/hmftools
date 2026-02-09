package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.expression.GeneCollectionSummary;

public class GcTranscriptCalculator
{
    private final IsofoxConfig mConfig;

    private final List<GeneData> mGeneDataList;

    private final Map<String,GcRatioCounts> mTranscriptGcRatioCache;

    private final GcRatioCounts mTranscriptFitGcCounts;
    private final double[] mGcRatioAdjustments;

    public GcTranscriptCalculator(final IsofoxConfig config)
    {
        mConfig = config;

        mGeneDataList = Lists.newArrayList();
        mTranscriptGcRatioCache = Maps.newHashMap();

        mTranscriptFitGcCounts = new GcRatioCounts();
        mGcRatioAdjustments = new double[mTranscriptFitGcCounts.size()];

        if(config.ExpGcRatiosFile != null)
            loadExpectedData();
    }

    public final double[] getGcRatioAdjustments() { return mGcRatioAdjustments; }
    public final GcRatioCounts getTranscriptFitGcCounts() { return mTranscriptFitGcCounts; }

    public void generateGcCountsFromFit(final List<GeneCollectionSummary> geneSummaries)
    {
        // use expected GC ratio counts and 1st-pass transcript fits to derive expected GC counts
        final double[] frequencies = mTranscriptFitGcCounts.getCounts();

        for(final GeneCollectionSummary geneSummary : geneSummaries)
        {
            final Map<String,Double> fitAllocations = geneSummary.getFitAllocations();

            for(Map.Entry<String,Double> entry : fitAllocations.entrySet())
            {
                final String transName = entry.getKey();

                double fitAlloc = entry.getValue();

                final GcRatioCounts transGcCounts = mTranscriptGcRatioCache.get(transName);

                if(transGcCounts == null)
                {
                    ISF_LOGGER.warn("genes({}) transcript({}) missing expected GC ratio counts from cache",
                            geneSummary.GeneNames, transName);
                    return;
                }

                final double[] transFrequencies = transGcCounts.getCounts();

                for(int i = 0; i < frequencies.length; ++i)
                {
                    frequencies[i] += transFrequencies[i] * fitAlloc;
                }
            }
        }
    }

    private static final double MAX_ADJUST_FACTOR = 3;
    private static final double MIN_ADJUST_FACTOR = 1 / MAX_ADJUST_FACTOR;
    private static final int ADJUST_LOWER_BOUND = 20;
    private static final int ADJUST_UPPER_BOUND = 80;

    public void calcGcRatioAdjustments(final GcRatioCounts globalGcCounts)
    {
        final double[] expectedFrequencies = mTranscriptFitGcCounts.getCounts();
        final double[] actualFrequencies = globalGcCounts.getCounts();

        double actualFrequencyTotal = globalGcCounts.getCountsTotal();
        double expectedFrequencyTotal = mTranscriptFitGcCounts.getCountsTotal();

        if(expectedFrequencyTotal == 0 || actualFrequencyTotal == 0)
        {
            ISF_LOGGER.error(format("invalid expected(%.0f) or actual(%.0f) totals", expectedFrequencyTotal, actualFrequencyTotal));
            return;
        }

        for(int i = 0; i < actualFrequencies.length; ++ i)
        {
            if(i < ADJUST_LOWER_BOUND || i > ADJUST_UPPER_BOUND)
            {
                mGcRatioAdjustments[i] = 1;
                continue;
            }

            double actualPerc = actualFrequencies[i] / actualFrequencyTotal;
            double expectedPerc = expectedFrequencies[i] / expectedFrequencyTotal;

            if(actualPerc > 0 && expectedPerc > 0)
                mGcRatioAdjustments[i] = max(min(expectedPerc/actualPerc, MAX_ADJUST_FACTOR), MIN_ADJUST_FACTOR);
            else if(actualPerc == 0 && expectedPerc == 0)
                mGcRatioAdjustments[i] = 1.0;
            else if(actualPerc == 0)
                mGcRatioAdjustments[i] = MAX_ADJUST_FACTOR;
            else
                mGcRatioAdjustments[i] = MIN_ADJUST_FACTOR;

            ISF_LOGGER.debug(format("ratio(%.2f) actual(%.6f) expected(%.6f) adjustment(%.3f)",
                    globalGcCounts.getRatios()[i], actualPerc, expectedPerc, mGcRatioAdjustments[i]));
        }
    }

    private void loadExpectedData()
    {
        if(!Files.exists(Paths.get(mConfig.ExpGcRatiosFile)))
        {
            ISF_LOGGER.error("invalid expected GC ratios file");
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.ExpGcRatiosFile));
            String fileDelim = inferFileDelimiter(mConfig.ExpGcRatiosFile);

            // skip field names
            String line = fileReader.readLine();

            if(line == null)
            {
                ISF_LOGGER.error("empty expected GC ratios file({})", mConfig.ExpGcRatiosFile);
                return;
            }

            GcRatioCounts gcRatioCounts = new GcRatioCounts();
            int expectedColCount = 1 + gcRatioCounts.size();

            while ((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim, -1);

                if(values.length != expectedColCount)
                {
                    ISF_LOGGER.error("invalid exp GC ratio data length({}) vs expected({}): {}", values.length, expectedColCount, line);
                    return;
                }

                final String transName = values[0];
                gcRatioCounts = new GcRatioCounts();
                final double[] frequencies = gcRatioCounts.getCounts();

                for(int i = 1; i < values.length; ++i)
                {
                    double counts = Double.parseDouble(values[i]);
                    frequencies[i - 1] = counts;
                }

                mTranscriptGcRatioCache.put(transName, gcRatioCounts);
            }

            ISF_LOGGER.info("loaded {} transcript expected GC ratios from file({})",
                    mTranscriptGcRatioCache.size(), mConfig.ExpGcRatiosFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load expected GC ratios file({}): {}", mConfig.ExpGcRatiosFile, e.toString());return;
        }
    }
}
