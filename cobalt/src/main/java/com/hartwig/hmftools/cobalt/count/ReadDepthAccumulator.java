package com.hartwig.hmftools.cobalt.count;

import static htsjdk.samtools.util.SequenceUtil.C;
import static htsjdk.samtools.util.SequenceUtil.G;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicIntegerArray;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;

// ReadDepthAccumulator accumulate read alignment blocks and calculate read depth per
// 1000 base windows.
//
// This class must be thread safe, maybe we should change it to non thread safe but put
// all updating of counts in one thread.
public class ReadDepthAccumulator
{
    // raw counts of read bases
    private static class ChromosomeWindowCounts
    {
        final String chromosome;
        final AtomicIntegerArray windowReadBaseCounts;
        final AtomicIntegerArray windowGcCounts;

        public ChromosomeWindowCounts(final String chromosome, final int numWindows)
        {
            this.chromosome = chromosome;
            this.windowReadBaseCounts = new AtomicIntegerArray(numWindows);
            this.windowGcCounts = new AtomicIntegerArray(numWindows);
        }

        public void addCount(int windowIndex, int count)
        {
            windowReadBaseCounts.getAndAdd(windowIndex, count);
        }
        public void addGcCount(int windowIndex, int count)
        {
            windowGcCounts.getAndAdd(windowIndex, count);
        }

        public int getCount(int windowIndex)
        {
            return windowReadBaseCounts.get(windowIndex);
        }
        public int getGcCount(int windowIndex)
        {
            return windowGcCounts.get(windowIndex);
        }
    }

    private final int mWindowSize;
    private final Map<String, ChromosomeWindowCounts> mChromosomeWindowCounts = new ConcurrentHashMap<>();

    public ReadDepthAccumulator(int windowSize)
    {
        mWindowSize = windowSize;
    }

    public void addChromosome(String chromosome, int chromosomeLength)
    {
        Validate.isTrue(!mChromosomeWindowCounts.containsKey(chromosome));
        int numWindows = chromosomeLength / mWindowSize;
        mChromosomeWindowCounts.put(chromosome, new ChromosomeWindowCounts(chromosome, numWindows));
    }

    @NotNull
    public List<DepthReading> getChromosomeReadDepths(String chromosome)
    {
        ChromosomeWindowCounts windowCounts = mChromosomeWindowCounts.get(chromosome);

        if(windowCounts == null)
        {
            // not a chromosome we keep track of
            return List.of();
        }

        List<DepthReading> readDepths = new ArrayList<>();

        for(int windowIndex = 0; windowIndex < windowCounts.windowReadBaseCounts.length(); ++windowIndex)
        {
            double basesCount = windowCounts.getCount(windowIndex);
            double depth = basesCount / mWindowSize;
            double gcContent = windowCounts.getGcCount(windowIndex) / basesCount;
            readDepths.add(new DepthReading(chromosome, getGenomePosition(windowIndex), depth, gcContent));
        }

        return readDepths;
    }

    // Add a read alignment to the base counts
    // this function is thread safe
    // genomeStart is 1 based and genomeEnd is inclusive
    // readStartIndex is 0 based
    public void addReadAlignmentToCounts(String chromosome, int genomeStart, int alignmentLength, byte[] readBases, int readStartIndex)
    {
        ChromosomeWindowCounts windowCounts = mChromosomeWindowCounts.get(chromosome);

        if(windowCounts == null)
        {
            // not a chromosome we keep track of
            return;
        }

        for(int windowIndex = getWindowIndex(genomeStart);; windowIndex++)
        {
            int windowStart = getGenomePosition(windowIndex);

            if(windowStart >= genomeStart + alignmentLength)
            {
                break;
            }

            if(windowIndex >= windowCounts.windowReadBaseCounts.length())
            {
                // this is possible as we omit the partial window at the end of chromosome
                break;
            }

            int numBasesInWindow = alignmentLength;

            Validate.isTrue(genomeStart < windowStart + mWindowSize);

            // if the start < window start, we have to account for it
            int startOffset = Math.max(0, windowStart - genomeStart);
            // if the window end is < alignment end
            int endOffset = Math.max(0, genomeStart + alignmentLength - windowStart - mWindowSize);

            numBasesInWindow -= startOffset;
            numBasesInWindow -= endOffset;

            Validate.isTrue(numBasesInWindow >= 0);
            windowCounts.addCount(windowIndex, numBasesInWindow);

            // next we need to count the GCs
            int numGcs = 0;
            for(int i = readStartIndex + startOffset; i < readStartIndex + startOffset + numBasesInWindow; ++i)
            {
                if(readBases[i] == G || readBases[i] == C)
                {
                    ++numGcs;
                }
            }

            windowCounts.addGcCount(windowIndex, numGcs);
        }
    }

    private int getWindowIndex(int position)
    {
        return (position - 1) / mWindowSize;
    }
    private int getGenomePosition(int windowIndex)
    {
        return windowIndex * mWindowSize + 1;
    }
}
