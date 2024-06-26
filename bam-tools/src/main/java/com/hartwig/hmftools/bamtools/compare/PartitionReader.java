package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.CompareUtils.basesMatch;
import static com.hartwig.hmftools.bamtools.compare.CompareUtils.stringsMatch;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.ORIG_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import static htsjdk.samtools.util.SequenceUtil.reverseComplement;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

public class PartitionReader implements Runnable
{
    private final String mName;
    private final CompareConfig mConfig;
    private final BamPartition mBamPartition;
    final BamReaderProvider mOrigBamReaderProvider;
    final BamReaderProvider mNewBamReaderProvider;

    private final Map<ReadKey, SAMRecord> mOrigBamReads = new HashMap<>();
    private final Map<ReadKey, SAMRecord> mNewBamReads = new HashMap<>();
    private final ReadWriter mReadWriter;

    @Nullable
    private final UnmatchedReadHandler mUnmatchedReadHandler;
    private final boolean mLogReadIds;

    private final Statistics mStats;

    public PartitionReader(
            final String name,
            final CompareConfig config, final BamPartition bamPartition,
            final BamReaderProvider bamAReaderProvider, final BamReaderProvider bamBReaderProvider,
            final ReadWriter readWriter, @Nullable final UnmatchedReadHandler unmatchedReadHandler)
    {
        mName = name;
        mConfig = config;
        mBamPartition = bamPartition;
        mOrigBamReaderProvider = bamAReaderProvider;
        mNewBamReaderProvider = bamBReaderProvider;
        mReadWriter = readWriter;
        mUnmatchedReadHandler = unmatchedReadHandler;
        mStats = new Statistics();
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public Statistics stats() { return mStats; }

    public void run()
    {
        BT_LOGGER.debug("processing {}", mName);

        // we process the records partition by partition
        // reads are stored inside a hash table and looked up by the read id
        try(final SAMRecordIterator origBamItr = mBamPartition.iterator(mOrigBamReaderProvider.getBamReader());
            final SAMRecordIterator newBamItr = mBamPartition.iterator(mNewBamReaderProvider.getBamReader()))
        {
            long origBamReadCount = 0;
            long newBamReadCount = 0;
            int origReadAlignStart = -1;
            int newReadAlignStart = -1;
            while(origBamItr.hasNext() || newBamItr.hasNext())
            {
                // decide which one to advance, or both
                boolean advanceOrig;
                boolean advanceNew;

                if(origReadAlignStart == newReadAlignStart && origBamItr.hasNext() && newBamItr.hasNext())
                {
                    // if the align start are the same try to advance both. This logic is required such that processing
                    // of unmapped region is not massively slowed down by only advancing one side all the time
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    // the one with smaller alignment start is advanced first
                    // this helps to process the same genomic location together
                    advanceOrig = !newBamItr.hasNext() || (origBamItr.hasNext() && (origReadAlignStart < newReadAlignStart));
                    advanceNew = !advanceOrig;
                }

                if(advanceOrig)
                {
                    origBamReadCount++;
                    final SAMRecord origBamRead = origBamItr.next();
                    processOrigBamRecord(origBamRead);
                    origReadAlignStart = origBamRead.getAlignmentStart();
                    // BT_LOGGER.info("advance orig");
                }
                if(advanceNew)
                {
                    newBamReadCount++;
                    final SAMRecord newBamRead = newBamItr.next();
                    processNewBamRecord(newBamRead);
                    newReadAlignStart = newBamRead.getAlignmentStart();
                    // BT_LOGGER.info("advance new");
                }

                // BT_LOGGER.printf(Level.INFO, "orig align start: %,d, new align start: %,d", origReadAlignStart, newReadAlignStart);

                // check if we need to dump the reads to unmatched read handler, this is to
                // protect against running out of memory
                if(mUnmatchedReadHandler != null &&
                        mOrigBamReads.size() + mNewBamReads.size() > mConfig.MaxCachedReadsPerThread)
                {
                    BT_LOGGER.trace("max cache read limit exceeded, writing reads to hash bam");
                    mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
                    mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
                    mOrigBamReads.clear();
                    mNewBamReads.clear();
                }
            }
            completePartition(origBamReadCount, newBamReadCount);
        }

        BT_LOGGER.printf(Level.DEBUG, "%s complete: origReads(%,d) newReads(%,d) diff(%,d)",
                mName, mStats.OrigReadCount, mStats.NewReadCount, mStats.DiffCount);
    }

    private void processOrigBamRecord(final SAMRecord origBamRead)
    {
        if(mLogReadIds && mConfig.LogReadIds.contains(origBamRead.getReadName()))
        {
            BT_LOGGER.debug("[orig bam] specific readId({})", origBamRead.getReadName());
        }

        if(excludeRead(origBamRead))
            return;

        /*
        if((mStats.RefReadCount % LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("partition({}) orig reads processed({})", mRegion, mStats.RefReadCount);
        }*/

        // create a key
        ReadKey readKey = ReadKey.from(origBamRead);

        // look it up see if can find corresponding beta record
        SAMRecord newBamRead = mNewBamReads.remove(readKey);

        if(newBamRead != null)
        {
            // found matching
            checkReadDetails(origBamRead, newBamRead);
        }
        else
        {
            // add to read map
            mOrigBamReads.put(readKey, origBamRead);
        }
    }

    private static final int LOG_COUNT = 1_000_000;

    private void processNewBamRecord(final SAMRecord newBamRead)
    {
        if(mLogReadIds && mConfig.LogReadIds.contains(newBamRead.getReadName()))
        {
            BT_LOGGER.debug("[new bam] specific readId({})", newBamRead.getReadName());
        }

        if(excludeRead(newBamRead))
            return;

        /*if((mStats.NewReadCount % LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("partition new reads processed({})", mStats.NewReadCount);
        }*/

        // create a key
        ReadKey readKey = ReadKey.from(newBamRead);

        // look it up see if can find corresponding beta record
        SAMRecord origRead = mOrigBamReads.remove(readKey);

        if(origRead != null)
        {
            // found matching
            checkReadDetails(origRead, newBamRead);
        }
        else
        {
            // add to beta read map
            mNewBamReads.put(readKey, newBamRead);
        }
    }

    private void completePartition(long origBamReads, long newBamReads)
    {
        /*BT_LOGGER.printf(Level.DEBUG, "partition(%s) complete, orig bam reads(%,d), new bam reads(%,d)",
                mName, origBamReads, newBamReads);*/

        if(mUnmatchedReadHandler != null)
        {
            // any reads not matched up is palmed off to the unmatched reads handler
            if(!mOrigBamReads.isEmpty())
            {
                mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
            }
            if(!mNewBamReads.isEmpty())
            {
                mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
            }
        }
        else
        {
            // write them out as mismatches
            mOrigBamReads.values().forEach(read -> mReadWriter.writeComparison(read, ORIG_ONLY, null));
            mStats.OrigReadCount += mOrigBamReads.size();
            mStats.DiffCount += mOrigBamReads.size();

            mNewBamReads.values().forEach(read -> mReadWriter.writeComparison(read, NEW_ONLY, null));
            mStats.NewReadCount += mNewBamReads.size();
            mStats.DiffCount += mNewBamReads.size();
        }

        mOrigBamReads.clear();
        mNewBamReads.clear();
    }

    private static final List<String> KEY_ATTRIBUTES = List.of(SUPPLEMENTARY_ATTRIBUTE, MATE_CIGAR_ATTRIBUTE);

    private void checkReadDetails(final SAMRecord read1, final SAMRecord read2)
    {
        ++mStats.OrigReadCount;
        ++mStats.NewReadCount;

        List<String> diffs = compareReads(read1, read2, mConfig);
        if(!diffs.isEmpty())
        {
            ++mStats.DiffCount;
            mReadWriter.writeComparison(read1, VALUE, diffs);
        }
    }

    static List<String> compareReads(final SAMRecord read1, final SAMRecord read2, final CompareConfig config)
    {
        List<String> diffs = new ArrayList<>();

        if(read1.getInferredInsertSize() != read2.getInferredInsertSize())
        {
            diffs.add(format("insertSize(%d/%d)", read1.getInferredInsertSize(), read2.getInferredInsertSize()));
        }

        if(read1.getMappingQuality() != read2.getMappingQuality())
        {
            diffs.add(format("mapQuality(%d/%d)", read1.getMappingQuality(), read2.getMappingQuality()));
        }

        if(!read1.getCigarString().equals(read2.getCigarString()))
        {
            diffs.add(format("cigar(%s/%s)", read1.getCigarString(), read2.getCigarString()));
        }

        if(read1.getFlags() != read2.getFlags())
        {
            if(read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag())
                diffs.add(format("negStrand(%s/%s)", read1.getReadNegativeStrandFlag(), read2.getReadNegativeStrandFlag()));

            if(!config.IgnoreDupDiffs && read1.getDuplicateReadFlag() != read2.getDuplicateReadFlag())
                diffs.add(format("duplicate(%s/%s)", read1.getDuplicateReadFlag(), read2.getDuplicateReadFlag()));
        }

        // check key attributes:
        for(String attribute : KEY_ATTRIBUTES)
        {
            if(attribute.equals(SUPPLEMENTARY_ATTRIBUTE) && config.IgnoreSupplementaryReads)
                continue;

            String readAttr1 = read1.getStringAttribute(attribute);
            String readAttr2 = read2.getStringAttribute(attribute);

            if(!Objects.equals(readAttr1, readAttr2)) // Objects.equals handles null case
            {
                diffs.add(format("attrib_%s(%s/%s)", attribute,
                        readAttr1 == null ? "missing" : readAttr1,
                        readAttr2 == null ? "missing" : readAttr2));
            }
        }

        // check the read bases, make sure we account for the read negative strand flag
        if(!basesMatch(read1.getReadString(), read1.getReadNegativeStrandFlag(),
                read2.getReadString(), read2.getReadNegativeStrandFlag()))
        {
            diffs.add(format("bases(%s/%s)",
                read1.getReadNegativeStrandFlag() ? reverseComplement(read1.getReadString()) : read1.getReadString(),
                read2.getReadNegativeStrandFlag() ? reverseComplement(read2.getReadString()) : read2.getReadString()));
        }
        // check the base qual, make sure we account for the read negative strand flag
        if(!stringsMatch(read1.getBaseQualityString(), read1.getReadNegativeStrandFlag(),
                read2.getBaseQualityString(), read2.getReadNegativeStrandFlag()))
        {
            diffs.add(format("baseQual(%s/%s)",
                read1.getReadNegativeStrandFlag() ? new StringBuilder(read1.getBaseQualityString()).reverse() : read1.getBaseQualityString(),
                read2.getReadNegativeStrandFlag() ? new StringBuilder(read2.getBaseQualityString()).reverse() : read2.getBaseQualityString()));
        }

        return diffs;
    }

    private boolean excludeRead(final SAMRecord read)
    {
        return read.isSecondaryAlignment() ||
              (mConfig.IgnoreAlterations && (read.hasAttribute(CONSENSUS_READ_ATTRIBUTE) || read.hasAttribute(UNMAP_ATTRIBUTE))) ||
              (mConfig.IgnoreConsensusReads && read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) ||
              (mConfig.IgnoreSupplementaryReads && read.getSupplementaryAlignmentFlag());
    }
}
