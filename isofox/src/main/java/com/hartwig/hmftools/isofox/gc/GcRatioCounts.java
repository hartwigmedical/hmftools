package com.hartwig.hmftools.isofox.gc;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GC_RATIO_BUCKET;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.isofox.IsofoxConfig;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class GcRatioCounts
{
    private final double[] mRatios; // the ratios
    private final double[] mCounts; // the frequencies of the ratios

    public static final double DEFAULT_GC_RATIO_BUCKET = 0.01;

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
        return ratio;
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
        double upperRatio = lowerRatio + GC_RATIO_BUCKET;

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
        if(otherCounts.length != mCounts.length)
            return;

        for(int i = 0; i < mCounts.length; ++i)
        {
            mCounts[i] += otherCounts[i];
        }
    }

    public double getPercentileRatio(double percentile)
    {
        double totalCounts = sumVector(mCounts);

        long currentTotal = 0;
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

    public static double calcGcRatioFromReadRegions(final IndexedFastaSequenceFile refFastaSeqFile, final String chromosome, final List<long[]> readRegions)
    {
        double gcRatioTotal = 0;
        int basesTotal = 0;
        for(final long[] region : readRegions)
        {
            final String bases = refFastaSeqFile.getSubsequenceAt(chromosome, region[SE_START], region[SE_END]).getBaseString();
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
            final BufferedWriter writer, final String name, final double[] data)
    {
        try
        {
            writer.write(name);

            for(int i = 0; i < data.length; ++i)
            {
                writer.write(String.format(",%.4f", data[i]));
            }

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write GC ratio file: {}", e.toString());
        }
    }

}
