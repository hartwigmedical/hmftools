package com.hartwig.hmftools.sig_analyser.loaders;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sig_analyser.common.CommonUtils.SIG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class PositionFrequencies
{
    private final Map<String, Map<Integer,Integer>> mChrPosBucketFrequencies;
    private final BufferedWriter mWriter;
    private final int mBucketSize;


    public PositionFrequencies(final String outputDir, final int bucketSize)
    {
        mBucketSize = bucketSize;
        mChrPosBucketFrequencies = Maps.newHashMap();
        mWriter = initialiseWriter(outputDir, mBucketSize);
    }

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
            SIG_LOGGER.error("failed to initialise position frequencies file: {}", e.toString());
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
            SIG_LOGGER.error("failed to write position frequency results: {}", e.toString());
        }
    }

}
