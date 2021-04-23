package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.sigs.VectorUtils.sumVectors;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GC_RATIO_BUCKET;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class GcRatioCounts
{
    private final double[] mRatios; // the ratios
    private final double[] mCounts; // the frequencies of the ratios

    public GcRatioCounts()
    {
        int itemCount = (int)floor(1 / GC_RATIO_BUCKET) + 1;
        mRatios = new double[itemCount];
        mCounts = new double[itemCount];

        double gcRatio = 0;

        for(int index = 0; index < mRatios.length; ++index)
        {
            mRatios[index] = gcRatio;
            gcRatio += GC_RATIO_BUCKET;
        }
    }

    public final double[] getRatios() { return mRatios; }
    public final double[] getCounts() { return mCounts; }
    public double getCountsTotal() { return sumVector(mCounts); }
    public int size() { return mRatios.length; }

    public static double roundGcRatio(double ratio)
    {
        return round(ratio/GC_RATIO_BUCKET) * GC_RATIO_BUCKET;
    }

    public static boolean isGC(char base) { return base == 'G' || base == 'C'; }

    public static int calcGcCount(final String bases)
    {
        int gcCount = 0;
        for (int i = 0; i < bases.length(); ++i)
        {
            if (isGC(bases.charAt(i)))
                ++gcCount;
        }

        return gcCount;
    }

    public static double calcGcRatio(final String bases)
    {
        int gcCount = calcGcCount(bases);
        return gcCount / (double) bases.length();
    }

    public void clearCounts()
    {
        for(int i = 0; i < mCounts.length; ++i)
        {
            mCounts[i] = 0;
        }
    }

    public void addGcRatio(double gcRatio)
    {
        final int[] gcRatioIndex = {-1, -1};
        final double[] counts = {1, 0};

        determineRatioData(gcRatio, gcRatioIndex, counts);

        addGcRatioCount(gcRatioIndex[0], counts[0]);

        if(gcRatioIndex[1] >= 0)
            addGcRatioCount(gcRatioIndex[1], counts[1]);
    }

    public void determineRatioData(double gcRatio, final int[] gcRatioIndex, final double[] counts)
    {
        double lowerRatio = floor(gcRatio/GC_RATIO_BUCKET) * GC_RATIO_BUCKET;

        if(lowerRatio >= 1)
        {
            gcRatioIndex[0] = getRatioIndex(lowerRatio);
            gcRatioIndex[1] = -1;
            counts[0] = 1;
            return;
        }

        // split proportionally amongst the 2 closest buckets
        gcRatioIndex[0] = getRatioIndex(lowerRatio);
        gcRatioIndex[1] = gcRatioIndex[0] + 1;
        counts[1] = (gcRatio - lowerRatio) / GC_RATIO_BUCKET;
        counts[0] = 1 - counts[1];
    }

    public void addGcRatioCount(int ratioIndex, double count)
    {
        if(ratioIndex < 0 || ratioIndex >= mCounts.length)
            return;

        mCounts[ratioIndex] += count;
    }

    private int getRatioIndex(double gcRatio)
    {
        double rawIndex = gcRatio * (mRatios.length - 1);
        return (int)min(max(0, round(rawIndex)), mRatios.length - 1);
    }

    public void mergeRatioCounts(final double[] otherCounts)
    {
        sumVectors(otherCounts, mCounts);
    }

    public double getPercentileRatio(double percentile)
    {
        double totalCounts = sumVector(mCounts);

        double currentTotal = 0;
        double prevRatio = 0;

        for(int i = 0; i < mRatios.length; ++i)
        {
            double gcRatio = mRatios[i];
            double frequency = mCounts[i];

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

    public static double calcGcRatioFromReadRegions(
            final RefGenomeInterface refGenome, final String chromosome, final List<int[]> readRegions)
    {
        double gcRatioTotal = 0;
        int basesTotal = 0;
        for(final int[] region : readRegions)
        {
            final String bases = refGenome.getBaseString(chromosome, region[SE_START], region[SE_END]);
            basesTotal += bases.length();
            gcRatioTotal += calcGcRatio(bases) * bases.length();
        }

        return gcRatioTotal / basesTotal;
    }

    public static BufferedWriter createReadGcRatioWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("gc_ratio_data.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneName");

            GcRatioCounts tmp = new GcRatioCounts();

            for(Double gcRatio : tmp.getRatios())
            {
                writer.write(String.format(",Gcr_%.2f", gcRatio));
            }

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
            final BufferedWriter writer, final String name, final double[] data, boolean asFractions)
    {
        try
        {
            writer.write(name);

            for(int i = 0; i < data.length; ++i)
            {
                if(asFractions)
                    writer.write(String.format(",%.4f", data[i]));
                else
                    writer.write(String.format(",%.0f", data[i]));
            }

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write GC ratio file: {}", e.toString());
        }
    }

}
