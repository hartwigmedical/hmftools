package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.UmiReadType.DUAL;
import static com.hartwig.hmftools.common.bam.UmiReadType.SINGLE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final String mUmiId; // the UMI if enabled

    private final FragmentCoords mFragmentCoords;
    private final List<SAMRecord> mReads;

    private SAMRecord mConsensusRead;
    private SAMRecord mPrimaryRead; // if no consensus is formed, the selected primary read
    private boolean mDualStrand;

    public DuplicateGroup(final String id, final SAMRecord read, final FragmentCoords fragmentCoords)
    {
        this(id, Lists.newArrayList(read), fragmentCoords);
    }

    public DuplicateGroup(final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        this(null, reads, fragmentCoords);
    }

    public DuplicateGroup(final String id, final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        mUmiId = id;
        mFragmentCoords = fragmentCoords;
        mReads = reads;

        mConsensusRead = null;
        mPrimaryRead = null;
        mDualStrand = false;
    }

    public void addRead(final SAMRecord read) { mReads.add(read); }
    public void addReads(final List<SAMRecord> reads) { mReads.addAll(reads); }

    public List<SAMRecord> reads() { return mReads; }
    public int readCount() { return mReads.size(); }

    public FragmentCoords fragmentCoordinates() { return mFragmentCoords; }

    public List<SAMRecord> duplicate() { return mReads; }
    public SAMRecord consensusRead() { return mConsensusRead; }

    public void setPrimaryRead(final SAMRecord read) { mPrimaryRead = read; }
    public boolean isPrimaryRead(final SAMRecord read) { return mPrimaryRead == read; }

    public String umiId() { return mUmiId; }

    public void registerDualStrand() { mDualStrand = true; }
    public boolean hasDualStrand() { return mDualStrand; }
    public String coordinatesKey() { return mFragmentCoords.Key; }

    public void formConsensusRead(final ConsensusReads consensusReads)
    {
        try
        {
            ConsensusReadInfo consensusReadInfo = consensusReads.createConsensusRead(mReads, mUmiId);

            // set consensus read attributes
            int firstInPairCount = (int)mReads.stream().filter(x -> x.getFirstOfPairFlag()).count();
            int readCount = mReads.size();
            boolean isPrimaryGroup = firstInPairCount >= readCount / 2;

            if(!isPrimaryGroup)
                firstInPairCount = readCount - firstInPairCount; // adjusted so both reads report the same ratio

            UmiReadType umiReadType = mDualStrand ? DUAL : SINGLE;

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount, firstInPairCount,  umiReadType);

            mConsensusRead = consensusReadInfo.ConsensusRead;
        }
        catch(Exception e)
        {
            RD_LOGGER.error("error forming consensus: {}", toString());

            for(SAMRecord read : mReads)
            {
                RD_LOGGER.error("read: {}", readToString(read));
            }

            e.printStackTrace();
        }
    }

    public String toString()
    {
        return format("id(%s) reads(%d) coords(%s)", mUmiId, mReads.size(), mFragmentCoords.Key);
    }
}
