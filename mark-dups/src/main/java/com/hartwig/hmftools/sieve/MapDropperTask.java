package com.hartwig.hmftools.sieve;

import static com.hartwig.hmftools.sieve.MapDropperConfig.MD_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.apache.commons.lang3.NotImplementedException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class MapDropperTask implements Callable
{
    private final MapDropperConfig mConfig;
    private final RecordWriter mWriter;
    private final SamReader mSamReader;
    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;
    private int mPrimaryUnmappedCounter;
    private int mSupplementaryDroppedCounter;

    public MapDropperTask(@NotNull final MapDropperConfig config, @NotNull RecordWriter writer)
    {
        mConfig = config;
        mWriter = writer;
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
        mPerfCounter = new PerformanceCounter("MapDropperTask");
        mReadCounter = 0;
        mPrimaryUnmappedCounter = 0;
        mSupplementaryDroppedCounter = 0;
    }

    @Override
    public Long call()
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            // TODO(m_cooper): Implement specific regions handling.
            MD_LOGGER.error("TODO: Implement specific regions handling.");
            throw new NotImplementedException("TODO: Implement specific regions handling.");
        }

        mPerfCounter.start();

        // Process reads in chunks which have the same read names.
        String lastReadName = null;
        List<SAMRecord> readNameRecords = null;
        final var it = mSamReader.iterator();
        while(it.hasNext())
        {
            final SAMRecord record = it.next();
            final String readName = record.getReadName();
            if(lastReadName == null)
            {
                lastReadName = readName;
                readNameRecords = new ArrayList<>();
                readNameRecords.add(record);
            }
            else if(lastReadName.equals(readName))
            {
                readNameRecords.add(record);
            }
            else
            {
                processReadGroup(readNameRecords);
                lastReadName = readName;
                readNameRecords = new ArrayList<>();
                readNameRecords.add(record);
            }
        }

        processReadGroup(readNameRecords);

        mPerfCounter.stop();
        MD_LOGGER.info("{} reads processed, {} primary reads unmapped, {} supplementary reads dropped", mReadCounter, mPrimaryUnmappedCounter, mSupplementaryDroppedCounter);
        return (long) 0;
    }

    private void unmapPrimaryRead(@NotNull final SAMRecord read, @NotNull final SAMRecord mate)
    {
        ++mPrimaryUnmappedCounter;

        read.setReadUnmappedFlag(true);
        read.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        read.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);

        mate.setMateUnmappedFlag(true);
        mate.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
    }

    private void processPrimaryAndSupplementaries(final boolean isFirstInPair, final @NotNull SAMRecord primaryRead,
            final @NotNull List<SAMRecord> supplementaryReads, final @NotNull SAMRecord primaryMate)
    {
        if(!primaryRead.getReadUnmappedFlag())
        {
            if(highFreqGMatch(primaryRead))
            {
                unmapPrimaryRead(primaryRead, primaryMate);
                mSupplementaryDroppedCounter += supplementaryReads.size();
                supplementaryReads.clear();
            }
            else
            {
                mSupplementaryDroppedCounter += supplementaryReads.stream().filter(this::highFreqGMatch).count();
                supplementaryReads.removeIf(this::highFreqGMatch);
            }
        }
        else if(!supplementaryReads.isEmpty())
        {
            MD_LOGGER.error("Primary read for {} in pair is unmapped, but it has supplementary reads.", isFirstInPair ? "first" : "second");
            System.exit(1);
        }
    }

    private void processReadGroup(@NotNull final List<SAMRecord> reads)
    {
        if(reads.isEmpty())
        {
            MD_LOGGER.error("Read group is empty.");
            System.exit(1);
        }

        mReadCounter += reads.size();
        final String readName = reads.get(0).getReadName();

        // Get first and second primary and supplementary reads.
        SAMRecord firstPrimary = null;
        SAMRecord secondPrimary = null;
        final List<SAMRecord> firstSupp = new ArrayList<>();
        final List<SAMRecord> secondSupp = new ArrayList<>();
        for(var record : reads)
        {
            if(record.getSupplementaryAlignmentFlag())
            {
                if(record.getFirstOfPairFlag())
                {
                    firstSupp.add(record);
                }
                else
                {
                    secondSupp.add(record);
                }
            }
            else
            {
                if(record.getFirstOfPairFlag())
                {
                    if(firstPrimary != null)
                    {
                        MD_LOGGER.error("Multiple primary reads for first in pair with read name {}", readName);
                        System.exit(1);
                    }

                    firstPrimary = record;
                }
                else
                {
                    if(secondPrimary != null)
                    {
                        MD_LOGGER.error("Multiple primary reads for second in pair with read name {}", readName);
                        System.exit(1);
                    }

                    secondPrimary = record;
                }
            }
        }

        if(firstPrimary == null)
        {
            MD_LOGGER.error("No primary read for first in pair with read name {}", readName);
            System.exit(1);
        }

        if(secondPrimary == null)
        {
            MD_LOGGER.error("No primary read for second in pair with read name {}", readName);
            System.exit(1);
        }

        // TODO(m_cooper): Handling past partition?
        processPrimaryAndSupplementaries(true, firstPrimary, firstSupp, secondPrimary);
        processPrimaryAndSupplementaries(false, secondPrimary, secondSupp, firstPrimary);

        mWriter.writeRead(firstPrimary);
        mWriter.writeRead(secondPrimary);
        mWriter.writeReads(firstSupp);
        mWriter.writeReads(secondSupp);
    }

    private boolean highFreqGMatch(@NotNull final SAMRecord read)
    {
        final String readBases = new String(read.getReadBases());
        final int baseStart = read.getAlignmentStart() - read.getUnclippedStart();
        final int baseEnd = read.getReadLength() - 1 - (read.getUnclippedEnd() - read.getAlignmentEnd());

        int gCount = 0;
        for(int i = baseStart; i <= baseEnd; ++i)
        {
            if(readBases.charAt(i) == 'G')
            {
                ++gCount;
            }
        }

        // TODO(m_cooper): Make fraction configurable.
        return gCount >= 0.9 * (baseEnd - baseStart + 1);
    }
}
