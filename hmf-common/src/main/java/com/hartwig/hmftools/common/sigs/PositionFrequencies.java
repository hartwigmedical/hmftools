package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.sigs.SigUtils.SU_LOGGER;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;

public class PositionFrequencies
{
    private final Map<String,Map<Integer,Integer>> mChrPosBucketFrequencies;
    private final BufferedWriter mWriter;
    private final int mBucketSize;

    public PositionFrequencies(final String outputDir, final int bucketSize)
    {
        mBucketSize = bucketSize;
        mChrPosBucketFrequencies = Maps.newHashMap();
        mWriter = initialiseWriter(outputDir, mBucketSize);
    }

    public final Map<String,Map<Integer,Integer>> getChrPosBucketFrequencies() { return mChrPosBucketFrequencies; }

    public void addPosition(final String chromosome, int position)
    {
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

    public void clear()
    {
        mChrPosBucketFrequencies.clear();
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }

    private static BufferedWriter initialiseWriter(final String outputDir, int bucketSize)
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

    public void writeResults(final String sampleId)
    {
        try
        {
            for(Map.Entry<String, Map<Integer,Integer>> chrEntry : mChrPosBucketFrequencies.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer,Integer> posEntry : chrEntry.getValue().entrySet())
                {
                    mWriter.write(String.format("%s,%s,%d,%d", sampleId, chromosome, posEntry.getKey(), posEntry.getValue()));
                    mWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            SU_LOGGER.error("failed to write position frequency results: {}", e.toString());
        }
    }

}
