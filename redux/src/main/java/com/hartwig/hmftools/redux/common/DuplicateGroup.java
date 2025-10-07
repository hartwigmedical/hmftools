package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNSET_COUNT;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.isIllumina;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.consensus.IlluminaRoutines.calculatePCRClusterCount;

import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final String mUmiId; // the UMI if enabled

    // with duplicate group collapsing some reads in mReads may not have mFragmentCoords FragmentCoords
    private final FragmentCoords mFragmentCoords;

    // contains reads that have potentially been merged due to UMI merging, all of these will be used for consensus read building
    private final List<SAMRecord> mReads;

    // contains reads that have been merged into this group due to jitter, these will not be used for consensus building
    private final List<SAMRecord> mNonConsensusReads;

    // contains reads that have been merged into this group due to poly-g UMI tail, these will not be used for consensus building
    private final List<SAMRecord> mPolyGUmiReads;

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
        mNonConsensusReads = Lists.newArrayList();
        mPolyGUmiReads = Lists.newArrayList();

        mConsensusRead = null;
        mPrimaryRead = null;
        mDualStrand = false;
    }

    public void addRead(final SAMRecord read) { mReads.add(read); }
    public void addReads(final List<SAMRecord> reads) { mReads.addAll(reads); }
    public void addNonConsensusReads(final List<SAMRecord> reads) { mNonConsensusReads.addAll(reads); }
    public void addPolyGUmiReads(final List<SAMRecord> reads) { mPolyGUmiReads.addAll(reads); }

    public List<SAMRecord> reads() { return mReads; }
    public List<SAMRecord> nonConsensusReads() { return mNonConsensusReads; }
    public List<SAMRecord> polyGUmiReads() { return mPolyGUmiReads; }

    public List<SAMRecord> allReads()
    {
        return Stream.concat(Stream.concat(mReads.stream(), mNonConsensusReads.stream()), mPolyGUmiReads.stream()).toList();
    }

    public int readCount() { return mReads.size() + mNonConsensusReads.size() + mPolyGUmiReads.size(); }

    public FragmentCoords fragmentCoordinates() { return mFragmentCoords; }

    public List<SAMRecord> duplicate() { return mReads; }
    public SAMRecord consensusRead() { return mConsensusRead; }

    public SAMRecord primaryRead() { return mPrimaryRead; }
    public void setPrimaryRead(final SAMRecord read) { mPrimaryRead = read; }
    public boolean isPrimaryRead(final SAMRecord read) { return mPrimaryRead == read; }

    public String umiId() { return mUmiId; }

    public void registerDualStrand() { mDualStrand = true; }
    public boolean hasDualStrand() { return mDualStrand; }

    public void formConsensusRead(final ConsensusReads consensusReads)
    {
        try
        {
            ConsensusReadInfo consensusReadInfo = consensusReads.createConsensusRead(mReads, mFragmentCoords, mUmiId);

            // set consensus read attributes
            int nonPolyGFirstInPairCount = (int) Stream.concat(mReads.stream(), mNonConsensusReads.stream())
                    .filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag())
                    .count();
            int nonPolyGReadCount = mReads.size() + mNonConsensusReads.size();
            boolean isPrimaryGroup = nonPolyGFirstInPairCount >= 0.5f * nonPolyGReadCount;

            if(!isPrimaryGroup)
                nonPolyGFirstInPairCount = nonPolyGReadCount - nonPolyGFirstInPairCount; // adjusted so both reads report the same ratio

            ConsensusType consensusType = mDualStrand ? DUAL : SINGLE;

            int pcrClusterCount = isIllumina() ? calculatePCRClusterCount(this) : UNSET_COUNT;

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount(), nonPolyGFirstInPairCount, consensusType, pcrClusterCount);

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
        return format("id(%s) reads(%d) coords(%s)", mUmiId, readCount(), mFragmentCoords.Key);
    }
}
