package com.hartwig.hmftools.bamtools.checker;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.READ_GROUP_ATTRIBUTE;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamReadLite;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private final String mReadId;
    private List<SAMRecord> mReads;

    private String mReadGroupIdId;
    private List<BamReadLite> mLiteReads;

    // counts of supplementary and secondary reads
    private short mExpectedSupplementaryCount;
    private short mReceivedSupplementaryCount;

    private String mFirstPrimaryCigar;
    private String mSecondPrimaryCigar;

    public Fragment(final SAMRecord read)
    {
        mReadId = read.getReadName();
        mReads = Lists.newArrayListWithExpectedSize(2);
        mLiteReads = null;
        mExpectedSupplementaryCount = 0;
        mReceivedSupplementaryCount = 0;
        mFirstPrimaryCigar = null;
        mSecondPrimaryCigar = null;

        addRead(read);
    }

    public String readId() { return mReadId; }
    public List<SAMRecord> reads() { return mReads != null ? mReads : Collections.emptyList(); }
    public int readCount() { return mReads != null ? mReads.size() : (mLiteReads != null ? mLiteReads.size() : 0); }
    public short expectedSupplementaryCount() { return mExpectedSupplementaryCount; }
    public short receivedSupplementaryCount() { return mReceivedSupplementaryCount; }

    public String firstPrimaryCigar() { return mFirstPrimaryCigar; }
    public String secondPrimaryCigar() { return mSecondPrimaryCigar; }

    public void addRead(final SAMRecord read)
    {
        if(mReads == null)
            mReads = Lists.newArrayListWithExpectedSize(2);

        mReads.add(read);

        if(read.getSupplementaryAlignmentFlag())
        {
            ++mReceivedSupplementaryCount;
        }
        else
        {
            if(read.getFirstOfPairFlag())
            {
                mFirstPrimaryCigar = read.getCigarString();
            }
            else
            {
                mSecondPrimaryCigar = read.getCigarString();
            }

            List<SupplementaryReadData> suppDataList = SupplementaryReadData.extractAlignments(read);

            if(suppDataList != null)
                mExpectedSupplementaryCount += suppDataList.size();
        }
    }

    public boolean hasPrimaryInfo() { return mFirstPrimaryCigar != null && mSecondPrimaryCigar != null; }

    public boolean isComplete()
    {
        return hasPrimaryInfo() && mExpectedSupplementaryCount == mReceivedSupplementaryCount;
    }

    public synchronized boolean mergeFragment(final Fragment fragment, final SAMFileHeader samFileHeader, final List<SAMRecord> completeReads)
    {
        // copy reads or primary cigar info if the fragment has already written the primaries
        transfer(fragment);

        // keep the fragment's reads if one or both primaries are missing
        if(!hasPrimaryInfo())
        {
            serialiseReads();
            return false;
        }

        // primary cigar info is complete, so write any cached reads
        deserialiseReads(samFileHeader);

        List<SAMRecord> fragCompleteReads = extractCompleteReads();

        completeReads.addAll(fragCompleteReads);

        // keep just its primary info if its is waiting on supplementaries only
        if(isComplete())
            return true;

        serialiseReads();
        return false;
    }

    public void transfer(final Fragment other)
    {
        other.reads().forEach(x -> addRead(x));

        if(other.hasPrimaryInfo() && !hasPrimaryInfo())
        {
            mFirstPrimaryCigar = other.firstPrimaryCigar();
            mSecondPrimaryCigar = other.secondPrimaryCigar();
            mExpectedSupplementaryCount += other.expectedSupplementaryCount();
        }
    }

    public List<SAMRecord> extractCompleteReads()
    {
        if(!hasPrimaryInfo() || mReads == null || mReads.isEmpty())
            return Collections.emptyList();

        setMateCigar();
        List<SAMRecord> reads = Lists.newArrayList(mReads);
        mReads.clear();
        return reads;
    }

    private void setMateCigar()
    {
        for(SAMRecord read : mReads)
        {
            if(read.getFirstOfPairFlag())
            {
                read.setAttribute(MATE_CIGAR_ATTRIBUTE, mSecondPrimaryCigar);
            }
            else
            {
                read.setAttribute(MATE_CIGAR_ATTRIBUTE, mFirstPrimaryCigar);
            }
        }
    }

    public String toString()
    {
        return format("id(%s) reads(%d) primary(first=%s second=%s) supps(expected=%d received=%d)",
                mReadId, readCount(), mFirstPrimaryCigar != null, mSecondPrimaryCigar != null,
                mExpectedSupplementaryCount, mReceivedSupplementaryCount);
    }

    public synchronized void serialiseReads()
    {
        if(mReads == null || mReads.isEmpty())
            return;

        if(mLiteReads == null)
            mLiteReads = Lists.newArrayListWithExpectedSize(mReads.size());

        if(mReadGroupIdId == null)
            mReadGroupIdId = mReads.get(0).getStringAttribute(READ_GROUP_ATTRIBUTE);

        for(SAMRecord read : mReads)
        {
            mLiteReads.add(new BamReadLite(read, true));
        }

        mReads.clear();
        mReads = null;
    }

    public void deserialiseReads(final SAMFileHeader samFileHeader)
    {
        if(mLiteReads == null || mLiteReads.isEmpty())
            return;

        if(mReads == null)
            mReads = Lists.newArrayListWithExpectedSize(mLiteReads.size());

        for(BamReadLite readLite : mLiteReads)
        {
            SAMRecord read = BamReadLite.from(readLite, samFileHeader, mReadId, mReadGroupIdId);
            mReads.add(read);
        }

        mLiteReads.clear();
        mLiteReads = null;
    }
}
