package com.hartwig.hmftools.common.sigs;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sigs.SigUtils.SU_LOGGER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;

public class PositionFrequencies
{
    private final Map<String,Map<Integer,Integer>> mChrPosBucketFrequencies;

    private final int mBucketSize;

    // position mapping
    private final Map<String,Integer> mChromosomeLengths;
    private final Map<String,Integer> mChromosomePosIndex;
    private int mPositionCacheSize;
    private int mMaxSampleCount;
    private final int[] mCounts;

    public static final int DEFAULT_POS_FREQ_BUCKET_SIZE = 500000;
    public static final int DEFAULT_POS_FREQ_MAX_SAMPLE_COUNT = 20000;

    public PositionFrequencies(final int bucketSize, final int maxSampleCount)
    {
        mBucketSize = bucketSize;
        mMaxSampleCount = maxSampleCount;

        mChromosomeLengths = Maps.newHashMap();
        mChromosomePosIndex = Maps.newHashMap();
        mChrPosBucketFrequencies = Maps.newHashMap();

        mPositionCacheSize = initialisePositionCache(mBucketSize, mChromosomeLengths, mChromosomePosIndex);

        mCounts = new int[mPositionCacheSize];
    }

    public int getBucketSize() { return mBucketSize; }
    public int getBucketCount() { return mPositionCacheSize; }
    public int getMaxSampleCount() { return mMaxSampleCount; }
    public final int[] getCounts() { return mCounts; }

    public final Map<String,Map<Integer,Integer>> getChrPosBucketFrequencies() { return mChrPosBucketFrequencies; }

    public void clear()
    {
        mChrPosBucketFrequencies.clear();

        for(int i = 0; i < mCounts.length; ++i)
        {
            mCounts[i] = 0;
        }
    }

    public boolean isValidChromosome(final String chromosome)
    {
        return mChromosomePosIndex.containsKey(chromosome);
    }

    public void addPosition(final String chromosome, int position)
    {
        int bucketIndex = getBucketIndex(mBucketSize, mChromosomePosIndex, chromosome, position);

        if(bucketIndex >= 0 && bucketIndex < mCounts.length)
            ++mCounts[bucketIndex];

        Map<Integer,Integer> positionMap = mChrPosBucketFrequencies.get(chromosome);

        if(positionMap == null)
        {
            positionMap = Maps.newHashMap();
            mChrPosBucketFrequencies.put(chromosome, positionMap);
        }

        int positionBucket = position / mBucketSize * mBucketSize;
        Integer frequency = positionMap.get(positionBucket);

        if(frequency == null)
            positionMap.put(positionBucket, 1);
        else
            positionMap.put(positionBucket, frequency + 1);
    }

    public static int initialisePositionCache(
            int bucketSize, final Map<String,Integer> chromosomeLengths, final Map<String,Integer> chrPosIndexMap)
    {
        int positionCacheSize = 0;

        if(bucketSize == 0)
            return positionCacheSize;

        chromosomeLengths.clear();

        final RefGenomeCoordinates refGenome37 = RefGenomeCoordinates.COORDS_37;
        final RefGenomeCoordinates refGenome38 = RefGenomeCoordinates.COORDS_38;

        int lastEndPosIndex = -1;

        for(HumanChromosome chr : HumanChromosome.values())
        {
            final String chromosome = chr.toString();
            int length = max(refGenome37.lengths().get(chr).intValue(), refGenome38.lengths().get(chr).intValue());
            chromosomeLengths.put(chromosome, length);

            // chromosomes will have position indices as: chr1 0-9, chr2 10-20 etc
            int startPosIndex = lastEndPosIndex > 0 ? lastEndPosIndex + 1 : 0;
            chrPosIndexMap.put(chromosome, startPosIndex);

            int positionCount = (int)ceil(length/(double)bucketSize);
            positionCacheSize += positionCount;

            lastEndPosIndex = startPosIndex + positionCount - 1;
        }

        return positionCacheSize;
    }

    public static final int INVALID_BUCKET_INDEX = -1;

    public static int getBucketIndex(
            final int bucketSize, final Map<String,Integer> chrPosIndexMap, final String chromosome, int position)
    {
        Integer chromosomePosIndex = chrPosIndexMap.get(chromosome);

        if(chromosomePosIndex == null)
            return INVALID_BUCKET_INDEX;

        int posBucket = (int)floor(position/(double)bucketSize);
        return chromosomePosIndex + posBucket;
    }

    public static String getChromosomeFromIndex(final Map<String,Integer> chrPosIndexMap, int bucketIndex)
    {
        int lastChrStartIndex = -1;
        String lastChromosome = "";
        for(HumanChromosome chr : HumanChromosome.values())
        {
            final String chromosome = chr.toString();
            int chrStartIndex = chrPosIndexMap.get(chromosome);

            if(lastChrStartIndex >= 0)
            {
                if(bucketIndex >= lastChrStartIndex && bucketIndex < chrStartIndex)
                {
                    return lastChromosome;
                }
            }

            lastChrStartIndex = chrStartIndex;
            lastChromosome = chromosome;
        }

        return lastChromosome;
    }

    public static int getPositionFromIndex(final Map<String,Integer> chrPosIndexMap, final String chromosome, int bucketIndex, int bucketSize)
    {
        int chrStartIndex = chrPosIndexMap.get(chromosome);
        return (bucketIndex - chrStartIndex) * bucketSize;
    }

    public static BufferedWriter createFrequencyCountsWriter(final String outputDir, int bucketSize)
    {
        if(outputDir == null || outputDir.isEmpty())
            return null;

        try
        {
            final String outputFile = String.format("%ssnv_position_freq_%d.csv", outputDir, bucketSize);
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("SampleId,Chromosome,Position,Count");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            SU_LOGGER.error("failed to initialise position frequencies file: {}", e.toString());
            return null;
        }
    }

    public void writeFrequencyCounts(final BufferedWriter writer, final String sampleId)
    {
        try
        {
            for(Map.Entry<String, Map<Integer,Integer>> chrEntry : mChrPosBucketFrequencies.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer,Integer> posEntry : chrEntry.getValue().entrySet())
                {
                    writer.write(String.format("%s,%s,%d,%d", sampleId, chromosome, posEntry.getKey(), posEntry.getValue()));
                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            SU_LOGGER.error("failed to write position frequency results: {}", e.toString());
        }
    }

}
