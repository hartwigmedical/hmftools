package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.CompareConfig.FULLY_UNMAPPED_PARTITION;
import static com.hartwig.hmftools.bamtools.compare.CompareUtils.compareReads;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.ORIG_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PartitionReader implements Callable<Void>
{
    private final String mName;
    private final CompareConfig mConfig;

    private final File mOrigBamFile;
    private final File mNewBamFile;
    private final ChrBaseRegion mPartition;

    private final Map<ReadKey,SAMRecord> mOrigBamReads;
    private final Map<ReadKey,SAMRecord> mNewBamReads;
    private final ReadWriter mReadWriter;

    @Nullable
    private final UnmatchedReadHandler mUnmatchedReadHandler;
    private final boolean mLogReadIds;

    @Nullable
    private final IgnoreRegionIndex mIgnoreRegionIndex;
    @Nullable
    private final Set<String> mIgnoredReadNames;

    private long mReadCount;
    private long mNextLogReadCount;

    private final Statistics mStats;

    public PartitionReader(
            final String name, final CompareConfig config, @Nullable final ChrBaseRegion bamPartition,
            final File origBamFile, final File newBamFile,
            final ReadWriter readWriter, @Nullable final UnmatchedReadHandler unmatchedReadHandler,
            @Nullable final IgnoreRegionIndex ignoreRegionIndex, @Nullable final Set<String> ignoredReadNames)
    {
        mName = name;
        mConfig = config;
        mPartition = bamPartition;

        mOrigBamFile = origBamFile;
        mNewBamFile = newBamFile;
        mReadWriter = readWriter;
        mUnmatchedReadHandler = unmatchedReadHandler;
        mIgnoreRegionIndex = ignoreRegionIndex;
        mIgnoredReadNames = ignoredReadNames;
        mStats = new Statistics();

        mOrigBamReads = new HashMap<>();
        mNewBamReads = new HashMap<>();

        mReadCount = 0;
        mNextLogReadCount = LOG_COUNT;

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public Statistics stats() { return mStats; }

    private String partitionStr() { return mPartition == null ? mName : mPartition.toString(); }

    @Override
    public Void call()
    {
        BT_LOGGER.trace("processing {}", partitionStr());

        SamReader origSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mOrigBamFile);

        SamReader newSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mNewBamFile);

        // reads are stored inside a hash table and looked up by the read ID
        SAMRecordIterator origBamIter;
        SAMRecordIterator newBamIter;

        boolean fullyUnmappedReads = mName.equals(FULLY_UNMAPPED_PARTITION);

        if(mPartition != null)
        {
            origBamIter = origSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
            newBamIter = newSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
        }
        else
        {
            if(fullyUnmappedReads)
            {
                origBamIter = origSamReader.queryUnmapped();
                newBamIter = newSamReader.queryUnmapped();
            }
            else
            {
                // all records (in the hash BAM)
                origBamIter = origSamReader.iterator();
                newBamIter = newSamReader.iterator();
            }
        }

        int origReadAlignStart = -1;
        int newReadAlignStart = -1;

        ReadKey origReadKey = null;
        SAMRecord origRead = null;
        ReadKey newReadKey = null;
        SAMRecord newRead = null;

        while(origBamIter.hasNext() || newBamIter.hasNext())
        {
            // decide which one to advance, or both
            boolean advanceOrig;
            boolean advanceNew;

            if(fullyUnmappedReads)
            {
                // while both BAMs have records, move ahead in one by a fixed amount to try to match off as many reads as possible
                // even though they aren't ordered
                if(origBamIter.hasNext() && newBamIter.hasNext())
                {
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    advanceOrig = origBamIter.hasNext();
                    advanceNew = !advanceOrig;
                }
            }
            else
            {
                if(origReadAlignStart == newReadAlignStart && origBamIter.hasNext() && newBamIter.hasNext())
                {
                    // if the align start are the same try to advance both. This logic is required such that processing
                    // of unmapped region is not massively slowed down by only advancing one side all the time
                    advanceOrig = true;
                    advanceNew = true;
                }
                else
                {
                    // the one with smaller alignment start is advanced first, so as to process the same genomic location together
                    advanceOrig = !newBamIter.hasNext() || (origBamIter.hasNext() && (origReadAlignStart < newReadAlignStart));
                    advanceNew = !advanceOrig;
                }
            }

            if(!advanceNew && !advanceOrig)
                break; // logical assert

            // in the event of the 2 reads matching, it will end up storing the original, only to retrieve again when the new is processed
            // but for 2 dissimilar BAMs that is probably unlikely
            // skip the total-record counters in the hash-bam re-processing pass to avoid double-counting
            final boolean countTotals = mUnmatchedReadHandler != null;

            if(advanceOrig)
            {
                if(countTotals)
                    ++mStats.OrigReadCount;
                origRead = origBamIter.next();
                origReadKey = checkAndFormReadKey(origRead);
                origReadAlignStart = origRead.getAlignmentStart();
            }

            if(advanceNew)
            {
                if(countTotals)
                    ++mStats.NewReadCount;
                newRead = newBamIter.next();
                newReadKey = checkAndFormReadKey(newRead);
                newReadAlignStart = newRead.getAlignmentStart();
            }

            if(origReadKey != null && newReadKey != null && origReadKey.equals(newReadKey))
            {
                checkReadDetails(origRead, newRead);
                origRead = null;
                newRead = null;
            }
            else
            {
                // rather than either wait for the next read to be a match, process any unmatched read now.
                // a null ReadKey means checkAndFormReadKey rejected the record (excludeRead, partition
                // out-of-range, etc) — drop it here rather than letting it slip into the unmatched cache
                // with a null map key, where it would either get compared with another null-keyed record
                // or be emitted as a spurious ORIG_ONLY / NEW_ONLY at partition end.
                if(origRead != null)
                {
                    if(origReadKey != null)
                        checkMatch(origRead, origReadKey, true);
                    origRead = null;
                    origReadKey = null;
                }

                if(newRead != null)
                {
                    if(newReadKey != null)
                        checkMatch(newRead, newReadKey, false);
                    newRead = null;
                    newReadKey = null;
                }
            }

            // check if we need to dump the reads to unmatched read handler, this is to protect against running out of memory
            int readCacheSize = mOrigBamReads.size() + mNewBamReads.size();

            if(mConfig.MaxCachedReadsPerThread > 0 && mUnmatchedReadHandler != null && readCacheSize > mConfig.MaxCachedReadsPerThread)
            {
                BT_LOGGER.debug("partition({}): cached reads({}) exceeds limit, writing reads to hash bam",
                        partitionStr(), readCacheSize);

                mUnmatchedReadHandler.handleOrigBamReads(mOrigBamReads.values());
                mUnmatchedReadHandler.handleNewBamReads(mNewBamReads.values());
                mOrigBamReads.clear();
                mNewBamReads.clear();
            }
        }

        if(origRead != null && origReadKey != null)
            mOrigBamReads.put(origReadKey, origRead);

        if(newRead != null && newReadKey != null)
            mNewBamReads.put(newReadKey, newRead);

        completePartition();

        if(mStats.OrigReadCount > 0 || mStats.NewReadCount > 0)
        {
            BT_LOGGER.trace("partition({}) complete: origReads({}) newReads({}) matched({}) diff({})",
                    partitionStr(), mStats.OrigReadCount, mStats.NewReadCount, mStats.Matched, mStats.DiffCount);
        }

        try
        {
            origSamReader.close();
            newSamReader.close();
        }
        catch(Exception e)
        {
            BT_LOGGER.error("failed to close BAMs: {}", e.toString());
        }

        return null;
    }

    private ReadKey checkAndFormReadKey(final SAMRecord read)
    {
        if(mPartition != null && !mPartition.containsPosition(read.getAlignmentStart()))
            return null;

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName()))
        {
            BT_LOGGER.debug("orig bam: specific readId({})", read.getReadName());
        }

        if(excludeRead(read))
            return null;

        checkLogReadCounts();

        // create a key
        return ReadKey.from(read);
    }

    private void checkMatch(final SAMRecord read, final ReadKey readKey, boolean isOriginal)
    {
        SAMRecord matchedRead = isOriginal ? mNewBamReads.remove(readKey) : mOrigBamReads.remove(readKey);

        if(matchedRead != null)
        {
            SAMRecord origRead = isOriginal ? read : matchedRead;
            SAMRecord newRead = !isOriginal ? read : matchedRead;
            checkReadDetails(origRead, newRead);
        }
        else
        {
            // add to read map
            if(isOriginal)
                mOrigBamReads.put(readKey, read);
            else
                mNewBamReads.put(readKey, read);
        }
    }

    private void completePartition()
    {
        /*BT_LOGGER.printf(Level.DEBUG, "partition(%s) complete, orig bam reads(%,d), new bam reads(%,d)",
                mName, origBamReads, newBamReads);*/

        if(mUnmatchedReadHandler != null)
        {
            // unmatched reads are then handled by the unmatched reads handler
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
            mOrigBamReads.values().forEach(read ->
            {
                DiffBucket bucket = DiffBucket.classify(ORIG_ONLY, read, null, null);
                mReadWriter.writeComparison(read, ORIG_ONLY, null, bucket);
                mStats.recordOrigOnly();
            });
            mNewBamReads.values().forEach(read ->
            {
                DiffBucket bucket = DiffBucket.classify(NEW_ONLY, read, null, null);
                mReadWriter.writeComparison(read, NEW_ONLY, null, bucket);
                mStats.recordNewOnly();
            });
        }

        mOrigBamReads.clear();
        mNewBamReads.clear();
    }

    private void checkReadDetails(final SAMRecord origRead, final SAMRecord newRead)
    {
        ++mStats.Matched;

        // note that the ReadKey of both reads match
        if(mConfig.IgnoreReduxAlterations && (origRead.hasAttribute(UNMAP_ATTRIBUTE) || newRead.hasAttribute(UNMAP_ATTRIBUTE)))
            return;

        List<String> diffs = compareReads(origRead, newRead, mConfig);

        if(!diffs.isEmpty())
        {
            DiffBucket bucket = DiffBucket.classify(VALUE, origRead, newRead, diffs);
            mStats.recordValue();
            mReadWriter.writeComparison(origRead, VALUE, diffs, bucket);
        }
    }

    private boolean excludeRead(final SAMRecord read)
    {
        if(mIgnoredReadNames != null && mIgnoredReadNames.contains(read.getReadName()))
        {
            ++mStats.IgnoredReadRecords;
            return true;
        }

        if(mConfig.IgnoreSecondaryReads && read.isSecondaryAlignment())
            return true;

        if(mConfig.IgnoreReduxAlterations)
        {
            if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE) || read.hasAttribute(UNMAP_ATTRIBUTE))
                return true;
        }

        if(mConfig.IgnoreConsensusReads && read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
            return true;

        if(mConfig.IgnoreSupplementaryReads && read.getSupplementaryAlignmentFlag())
            return true;

        return false;
    }

    public void scanIgnoredReadNames()
    {
        if(mIgnoreRegionIndex == null || mIgnoredReadNames == null)
            return;

        SamReader origSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mOrigBamFile);

        SamReader newSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(mNewBamFile);

        SAMRecordIterator origBamIter;
        SAMRecordIterator newBamIter;

        boolean fullyUnmappedReads = mName.equals(FULLY_UNMAPPED_PARTITION);

        if(mPartition != null)
        {
            origBamIter = origSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
            newBamIter = newSamReader.query(mPartition.chromosome(), mPartition.start(), mPartition.end(), false);
        }
        else if(fullyUnmappedReads)
        {
            try { origSamReader.close(); newSamReader.close(); }
            catch(Exception e) { BT_LOGGER.error("failed to close BAMs"); }
            return;
        }
        else
        {
            origBamIter = origSamReader.iterator();
            newBamIter = newSamReader.iterator();
        }

        scanIterator(origBamIter);
        scanIterator(newBamIter);

        try
        {
            origSamReader.close();
            newSamReader.close();
        }
        catch(Exception e)
        {
            BT_LOGGER.error("failed to close BAMs: {}", e.toString());
        }
    }

    private void scanIterator(final SAMRecordIterator iter)
    {
        String lastFlaggedReadName = null;
        while(iter.hasNext())
        {
            final SAMRecord read = iter.next();

            if(read.getReadUnmappedFlag())
                continue;

            if(mPartition != null && !mPartition.containsPosition(read.getAlignmentStart()))
                continue;

            final String readName = read.getReadName();
            if(readName.equals(lastFlaggedReadName))
                continue;

            if(mIgnoreRegionIndex.overlapsAboveHalf(
                    read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd()))
            {
                mIgnoredReadNames.add(readName);
                lastFlaggedReadName = readName;
            }
        }
    }

    private static final int LOG_COUNT = 1_000_000;

    private void checkLogReadCounts()
    {
        ++mReadCount;

        if(mReadCount >= mNextLogReadCount)
        {
            mNextLogReadCount += LOG_COUNT;

            BT_LOGGER.debug("partition({}) reads processed total({}) orig({}) new({}) matched({}) diffs({}) cache(old={} new={})",
                    partitionStr(), mReadCount, mStats.OrigReadCount, mStats.NewReadCount, mStats.Matched, mStats.DiffCount,
                    mOrigBamReads.size(), mNewBamReads.size());
        }
    }
}
