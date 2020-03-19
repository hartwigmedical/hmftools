package com.hartwig.hmftools.isofox.gc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GC_RATIO_BUCKET;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class GcRatioCounts
{
    private final Map<Double,Integer> mRatioCounts; // frequency of each GC-ratio
    private final Map<Double,Integer> mGeneRatioCounts;

    public static final double DEFAULT_GC_RATIO_BUCKET = 0.01;

    public GcRatioCounts()
    {
        mRatioCounts = Maps.newHashMap();
        mGeneRatioCounts = Maps.newHashMap();
    }

    public void processRead(final String readBases)
    {
        double gcRatio = calcGcRatio(readBases);
        addGcRatioCount(gcRatio);
    }

    public static double roundGcRatio(double ratio)
    {
        return round(ratio/GC_RATIO_BUCKET) * GC_RATIO_BUCKET;
    }

    public static double calcGcRatio(final String bases)
    {
        int gcCount = 0;
        for (int i = 0; i < bases.length(); ++i)
        {
            if (bases.charAt(i) == 'C' || bases.charAt(i) == 'G')
                ++gcCount;
        }

        double ratio = gcCount / (double) bases.length();
        return roundGcRatio(ratio);
    }

    public final Map<Double,Integer> getRatioCounts() { return mRatioCounts; }
    public final Map<Double,Integer> getGeneRatioCounts() { return mGeneRatioCounts; }

    public void clearGeneCounts() { mGeneRatioCounts.clear(); }

    private void addGcRatioCount(double gcRatio)
    {
        addGcRatioCount(mRatioCounts, gcRatio);
        addGcRatioCount(mGeneRatioCounts, gcRatio);
    }

    private static void addGcRatioCount(final Map<Double,Integer> ratioCounts, double gcRatio)
    {
        addGcRatioCount(ratioCounts, gcRatio, 1);
    }

    private static void addGcRatioCount(final Map<Double,Integer> ratioCounts, double gcRatio, int count)
    {
        Integer ratioCount = ratioCounts.get(gcRatio);

        if (ratioCount == null)
        {
            ratioCounts.put(gcRatio, count);
        }
        else
        {
            ratioCounts.put(gcRatio, ratioCount + count);
        }
    }

    public void mergeRatioCounts(final Map<Double,Integer> otherCounts)
    {
        otherCounts.entrySet().forEach(x -> addGcRatioCount(mRatioCounts, x.getKey(), x.getValue()));
    }

    public double getPercentileRatio(double percentile)
    {
        double totalCounts = mRatioCounts.values().stream().mapToLong(x -> x.longValue()).sum();

        final List<Double> sortedRatios = Lists.newArrayList();

        for(Double gcRatio : mRatioCounts.keySet())
        {
            int index = 0;
            while(index < sortedRatios.size())
            {
                if(gcRatio < sortedRatios.get(index))
                    break;

                ++index;
            }

            sortedRatios.add(index, gcRatio);
        }

        long currentTotal = 0;
        double prevRatio = 0;

        for(Double gcRatio : sortedRatios)
        {
            int frequency = mRatioCounts.get(gcRatio);

            double nextPercTotal = (currentTotal + frequency) / totalCounts;

            if(nextPercTotal >= percentile)
            {
                double medianRatio = prevRatio > 0 ? (prevRatio + gcRatio) * 0.5 : gcRatio;
                return medianRatio;
            }

            currentTotal += frequency;
            prevRatio = gcRatio;

        }

        return 0;
    }

    public static BufferedWriter createReadGcRatioWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("read_gc_ratios.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,GcRatio,Count");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create GC ratio writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeReadGcRatioCounts(
            final BufferedWriter writer, final EnsemblGeneData geneData, final Map<Double,Integer> ratioCounts)
    {
        try
        {
            for(Map.Entry<Double,Integer> entry : ratioCounts.entrySet())
            {
                writer.write(String.format("%s,%s,%.4f,%d",
                        geneData != null ? geneData.GeneId : "ALL", geneData != null ? geneData.GeneName : "ALL",
                        entry.getKey(), entry.getValue()));

                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write GC ratio file: {}", e.toString());
        }
    }

}
