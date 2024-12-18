package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.redux.PartitionReader.fullyUnmapped;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class ReadChecker
{
    private final boolean mUseReadCache;

    private final Map<String,Map<String,ReadInfo>> mChrUnmappedReadMap;

    private final AtomicLong mReadRemapped;
    private final AtomicLong mReadFullyUnmapped;

    private final AtomicLong mCheckReadRemapped;
    private final AtomicLong mCheckReadFullyUnmapped;

    public ReadChecker(boolean cacheReads)
    {
        mUseReadCache = cacheReads;
        mChrUnmappedReadMap = Maps.newHashMap();

        mReadRemapped = new AtomicLong(0);
        mReadFullyUnmapped = new AtomicLong(0);
        mCheckReadRemapped = new AtomicLong(0);
        mCheckReadFullyUnmapped = new AtomicLong(0);
    }

    private enum ReadState
    {
        UNCHANGED,
        REMAPPED,
        FULLY_UNMAPPED;
    }

    private class ReadInfo
    {
        public final String Info;
        public final ReadState State;

        public ReadInfo(final String info, final ReadState state)
        {
            Info = info;
            State = state;
        }

        public String toString() { return format("%s: %s", State, Info); }
    }

    public void addUnmappedRead(final SAMRecord read, final String chromosome, final int readStart)
    {
        if(read.isSecondaryOrSupplementary())
            return;

        boolean fullyUnmapped = fullyUnmapped(read);

        if(fullyUnmapped)
            mReadFullyUnmapped.incrementAndGet();
        else
            mReadRemapped.incrementAndGet();

        if(!mUseReadCache)
            return;

        cacheUnmappedRead(read, chromosome, readStart, fullyUnmapped);
    }

    private synchronized void cacheUnmappedRead(final SAMRecord read, final String chromosome, int readStart, boolean fullyUnmapped)
    {
        String readInfo = format("%s:%d mate(%s:%d) %s",
                chromosome, readStart,
                read.getReadPairedFlag() ? read.getMateReferenceName() : NO_CHROMOSOME_NAME,
                read.getReadPairedFlag() ? read.getMateAlignmentStart() : NO_POSITION, read.getCigarString());

        Map<String,ReadInfo> chrEntries = mChrUnmappedReadMap.get(chromosome);

        if(chrEntries == null)
        {
            chrEntries = Maps.newHashMap();
            mChrUnmappedReadMap.put(chromosome, chrEntries);
        }

        String readId = formReadId(read);

        ReadState readState = fullyUnmapped ? ReadState.FULLY_UNMAPPED : ReadState.REMAPPED;
        chrEntries.put(readId, new ReadInfo(readInfo, readState));
    }

    private static String formReadId(final SAMRecord read)
    {
        return format("%s_%d", read.getReadName(), !read.getReadPairedFlag() || read.getFirstOfPairFlag() ? 1 : 2);
    }

    public void checkRead(final SAMRecord read, final String chromosome, int readStart, boolean fullyUnmapped)
    {
        if(read.isSecondaryOrSupplementary())
            return;

        if(read.getReadUnmappedFlag())
        {
            if(fullyUnmapped)
                mCheckReadFullyUnmapped.incrementAndGet();
            else
                mCheckReadRemapped.incrementAndGet();
        }

        if(!mUseReadCache)
            return;

        checkCachedRead(read, chromosome, readStart, fullyUnmapped);
    }

    private synchronized void checkCachedRead(final SAMRecord read, final String chromosome, int readStart, boolean fullyUnmapped)
    {
        boolean expectCached = read.getReadUnmappedFlag();
        String readId = formReadId(read);

        Map<String,ReadInfo> chrEntries = mChrUnmappedReadMap.get(chromosome);

        if(chrEntries == null)
        {
            RD_LOGGER.warn("read chr({}:{}) id({}) not in cache, details: {}",
                    chromosome, readStart, readId, readToString(read));

            return;
        }

        ReadInfo readInfo = chrEntries.remove(readId);

        if(expectCached)
        {
            if(readInfo == null)
            {
                RD_LOGGER.warn("read chr({}:{}) id({}) not in cache, details: {}",
                        chromosome, readStart, readId, readToString(read));
                return;
            }
        }
        else
        {
            if(readInfo != null)
            {
                RD_LOGGER.warn("read chr({}:{}) id({}) in cache when expected not, details: {}",
                        chromosome, readStart, readId, readToString(read));
            }

            return;
        }

        ReadState expectedReadState = !read.getReadUnmappedFlag() ? ReadState.UNCHANGED :
                (fullyUnmapped ? ReadState.FULLY_UNMAPPED : ReadState.REMAPPED);

        if(expectedReadState != readInfo.State)
        {
            RD_LOGGER.warn("read({}) state differs cached({}) new({})", readInfo, readInfo.State, expectedReadState);
        }
    }

    public void logUnmatchedUnmappedReads()
    {
        if(mReadRemapped.get() != mCheckReadRemapped.get() || mReadFullyUnmapped.get() != mCheckReadFullyUnmapped.get())
        {
            RD_LOGGER.warn("unmappped mismatch: unmapping(unmapped={} full={}) partition reads(unmapped={} full={})",
                    mReadRemapped.get(), mReadFullyUnmapped.get(), mCheckReadRemapped.get(), mCheckReadFullyUnmapped.get());
        }

        if(!mUseReadCache)
            return;

        for(Map<String,ReadInfo> chrMap : mChrUnmappedReadMap.values())
        {
            for(Map.Entry<String,ReadInfo> entry : chrMap.entrySet())
            {
                RD_LOGGER.warn("read({}) details({}) remains cached", entry.getKey(), entry.getValue());
            }
        }

        mChrUnmappedReadMap.clear();
    }
}
