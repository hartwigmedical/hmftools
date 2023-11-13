package com.hartwig.hmftools.cobalt.count;

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

        public ChromosomeWindowCounts(final String chromosome, final int numWindows)
        {
            this.chromosome = chromosome;
            this.windowReadBaseCounts = new AtomicIntegerArray(numWindows);
        }

        public void addCount(int windowIndex, int count)
        {
            windowReadBaseCounts.getAndAdd(windowIndex, count);
        }

        public int getCount(int windowIndex)
        {
            return windowReadBaseCounts.get(windowIndex);
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

        if(chromosomeLength % mWindowSize > 0.8 * mWindowSize)
        {
            // we add extra window if it is at least 0.8 * window size
            ++numWindows;
        }

        mChromosomeWindowCounts.put(chromosome, new ChromosomeWindowCounts(chromosome, numWindows));
    }

    @NotNull
    public List<ReadDepth> getChromosomeReadDepths(String chromosome)
    {
        ChromosomeWindowCounts windowCounts = mChromosomeWindowCounts.get(chromosome);

        if(windowCounts == null)
        {
            // not a chromosome we keep track of
            return List.of();
        }

        List<ReadDepth> readDepths = new ArrayList<>();

        for(int windowIndex = 0; windowIndex < windowCounts.windowReadBaseCounts.length(); ++windowIndex)
        {
            double depth = (double)windowCounts.getCount(windowIndex) / mWindowSize;
            readDepths.add(new ReadDepth(chromosome, getGenomePosition(windowIndex), depth));
        }

        return readDepths;
    }

    // Add a read alignment to the base counts
    // this function is thread safe
    public void addReadAlignmentToCounts(String chromosome, int start, int end)
    {
        ChromosomeWindowCounts windowCounts = mChromosomeWindowCounts.get(chromosome);

        if(windowCounts == null)
        {
            // not a chromosome we keep track of
            return;
        }

        for(int windowIndex = getWindowIndex(start);; windowIndex++)
        {
            int windowStart = getGenomePosition(windowIndex);

            if(windowStart > end)
            {
                break;
            }

            if(windowIndex >= windowCounts.windowReadBaseCounts.length())
            {
                break;
            }

            int countToAdd = mWindowSize;

            Validate.isTrue(start < windowStart + mWindowSize);

            // if the start does not cover window start, we have to account for it
            countToAdd -= Math.max(0, start - windowStart);

            // if the window end is > alignment end
            countToAdd -= Math.max(0, (windowStart + mWindowSize - 1) - end);

            Validate.isTrue(countToAdd >= 0);

            windowCounts.addCount(windowIndex, countToAdd);
        }
    }

    // get the raw base count of the window, useful for testing
    public int getWindowRawBaseCount(String chromosome, int genomePosition)
    {
        int windowIndex = getWindowIndex(genomePosition);

        ChromosomeWindowCounts windowCounts = mChromosomeWindowCounts.get(chromosome);

        if(windowCounts == null)
        {
            // not a chromosome we keep track of
            return 0;
        }

        if(windowIndex >= windowCounts.windowReadBaseCounts.length())
        {
            // over the end
            return 0;
        }

        return windowCounts.getCount(windowIndex);
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
