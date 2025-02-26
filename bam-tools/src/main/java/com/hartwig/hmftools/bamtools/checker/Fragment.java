package com.hartwig.hmftools.bamtools.checker;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private final String mReadId;
    private final List<SAMRecord> mReads;

    // counts of supplementary and secondary reads
    private int mExpectedSupplementaryCount;
    private int mReceivedSupplementaryCount;

    private String mFirstPrimaryCigar;
    private String mSecondPrimaryCigar;

    public Fragment(final SAMRecord read)
    {
        mReadId = read.getReadName();
        mReads = Lists.newArrayListWithExpectedSize(2);
        mExpectedSupplementaryCount = 0;
        mReceivedSupplementaryCount = 0;

        addRead(read);
    }

    public String readId() { return mReadId; }
    public List<SAMRecord> reads() { return mReads; }

    public void addRead(final SAMRecord read)
    {
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

    public List<SAMRecord> extractCompleteReads()
    {
        if(!hasPrimaryInfo() || mReads.isEmpty())
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
        return format("id(%s) primary(first=%s second=%s) supps(expected=%d received=%d)",
                mReadId, mFirstPrimaryCigar != null, mSecondPrimaryCigar != null, mExpectedSupplementaryCount, mReceivedSupplementaryCount);
    }
}
