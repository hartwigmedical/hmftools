package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.genome.window.Window;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// read bam file to find the count for a given chromsome
// always return a value for each window, even for window without any read
// it will still create a ReadCount of 0
public class ChromosomeReadCount implements Callable<ChromosomeReadCount>
{
    private final File mInputFile;
    private final SamReaderFactory mReaderFactory;
    private final Chromosome mChromosome;
    private final List<ReadCount> mResults;
    private final int mMinMappingQuality;
    private final Window mWindow;

    private int mStart;
    private int mCount;

    public ChromosomeReadCount(
            final File inputFile, final SamReaderFactory readerFactory, final Chromosome chromosome,
            final int windowSize, final int minMappingQuality)
    {
        mInputFile = inputFile;
        mReaderFactory = readerFactory;
        mChromosome = chromosome;
        mMinMappingQuality = minMappingQuality;
        mWindow = new Window(windowSize);
        mResults = Lists.newArrayList();

        mStart = 1;
        mCount = 0;
    }

    @Override
    public ChromosomeReadCount call() throws Exception
    {
        CB_LOGGER.info("Generating windows on chromosome {}", mChromosome);

        try(final SamReader reader = mReaderFactory.open(mInputFile))
        {
            final SAMRecordIterator iterator = reader.query(mChromosome.contig, 0, 0, true);
            while(iterator.hasNext())
            {
                addRecord(iterator.next());
            }
        }
        // finalise the last window
        addReadCount(mStart, mCount);
        return this;
    }

    public Chromosome chromosome()
    {
        return mChromosome;
    }

    public List<ReadCount> readCount()
    {
        return mResults;
    }

    private void addRecord(final SAMRecord record)
    {
        if(!isEligible(record))
            return;

        int window = windowPosition(record.getAlignmentStart());
        while (mStart != window)
        {
            addReadCount(mStart, mCount);
            mStart += mWindow.getSize();
            mCount = 0;
        }

        mCount++;
    }

    private void addReadCount(int position, int count)
    {
        mResults.add(ImmutableReadCount.builder().chromosome(mChromosome.contig).position(position).readCount(count).build());
    }

    private boolean isEligible(final SAMRecord record)
    {
        return record.getMappingQuality() >= mMinMappingQuality
                && !(record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.isSecondaryOrSupplementary());
    }

    private int windowPosition(int position)
    {
        return mWindow.start(position);
    }
}
