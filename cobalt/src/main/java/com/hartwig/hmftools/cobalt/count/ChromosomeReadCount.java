package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.window.Window;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ChromosomeReadCount implements Callable<ChromosomeReadCount>
{
    private final File mInputFile;
    private final SamReaderFactory mReaderFactory;
    private final String mChromosome;
    private final int mChromosomeLength;
    private final List<ReadCount> mResults;
    private final int mMinMappingQuality;
    private final Window mWindow;

    private int mStart;
    private int mCount;

    public ChromosomeReadCount(
            final File inputFile, final SamReaderFactory readerFactory, final String chromosome,
            final int chromosomeLength, final int windowSize, final int minMappingQuality)
    {
        mInputFile = inputFile;
        mReaderFactory = readerFactory;
        mChromosome = chromosome;
        mChromosomeLength = chromosomeLength;
        mMinMappingQuality = minMappingQuality;
        mWindow = new Window(windowSize);
        mResults = Lists.newArrayList();

        mStart = 1;
        mCount = -1;
    }

    @Override
    public ChromosomeReadCount call() throws Exception
    {
        CB_LOGGER.info("Generating windows on chromosome {}", mChromosome);

        try(final SamReader reader = mReaderFactory.open(mInputFile))
        {
            final SAMRecordIterator iterator = reader.query(mChromosome, 0, 0, true);
            while(iterator.hasNext())
            {
                addRecord(iterator.next());
            }
        }
        return this;
    }

    public Chromosome chromosome()
    {
        return HumanChromosome.fromString(mChromosome);
    }

    public List<ReadCount> readCount()
    {
        addReadCount(mStart, mCount);

        int lastWindowPosition = lastWindowPosition();
        if(mResults.get(mResults.size() - 1).position() < lastWindowPosition)
        {
            addReadCount(lastWindowPosition, -1);
        }

        return mResults;
    }

    private void addRecord(final SAMRecord record)
    {
        if(!isEligible(record))
            return;

        int window = windowPosition(record.getAlignmentStart());
        if(mStart != window)
        {
            addReadCount(mStart, mCount);
            mStart = window;
            mCount = 0;
        }

        mCount++;
    }

    private void addReadCount(int position, int count)
    {
        mResults.add(ImmutableReadCount.builder().chromosome(mChromosome).position(position).readCount(count).build());
    }

    private boolean isEligible(final SAMRecord record)
    {
        return record.getMappingQuality() >= mMinMappingQuality
                && !(record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.isSecondaryOrSupplementary());
    }

    private int lastWindowPosition()
    {
        return windowPosition(mChromosomeLength);
    }

    private int windowPosition(int position)
    {
        return mWindow.start(position);
    }
}
