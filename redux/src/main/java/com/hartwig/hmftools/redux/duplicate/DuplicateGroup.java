package com.hartwig.hmftools.redux.duplicate;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNSET_COUNT;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.isIllumina;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.consensus.IlluminaRoutines.calculatePCRClusterCount;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final String mUmi; // the UMI if enabled

    // with duplicate group collapsing some reads in mReads may not have mFragmentCoords FragmentCoords
    private final FragmentCoords mFragmentCoords;

    // contains raw duplicate reads and those merged due to UMI matching and these will be used for consensus
    private final List<SAMRecord> mReads;

    // contains reads that have been merged due to jitter (poly-G, UMI) but will not be used for consensus
    private final List<SAMRecord> mNonConsensusReads;

    // contains reads that have been merged into this group due to poly-g UMI tail, these will not be used for consensus building
    private final List<SAMRecord> mPolyGUmiReads;

    private SAMRecord mConsensusRead;
    private SAMRecord mPrimaryRead; // if no consensus is formed, the selected primary read
    private boolean mDualStrand;

    public static final Comparator<DuplicateGroup> DUPLICATE_GROUP_COMPARATOR = Comparator.comparingInt(DuplicateGroup::totalReadCount).reversed();

    public DuplicateGroup(final String id, final SAMRecord read, final FragmentCoords fragmentCoords)
    {
        this(id, Lists.newArrayList(read), fragmentCoords);
    }

    public DuplicateGroup(final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        this(null, reads, fragmentCoords);
    }

    public DuplicateGroup(final String umi, final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        mUmi = umi;
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

    public int totalReadCount() { return mReads.size() + mNonConsensusReads.size() + mPolyGUmiReads.size(); }

    public FragmentCoords fragmentCoordinates() { return mFragmentCoords; }

    public List<SAMRecord> duplicate() { return mReads; }
    public SAMRecord consensusRead() { return mConsensusRead; }

    public SAMRecord primaryRead() { return mPrimaryRead; }
    public void setPrimaryRead(final SAMRecord read) { mPrimaryRead = read; }
    public boolean isPrimaryRead(final SAMRecord read) { return mPrimaryRead == read; }

    public String umi() { return mUmi; }

    public void registerDualStrand() { mDualStrand = true; }
    public boolean hasDualStrand() { return mDualStrand; }

    public void formConsensusRead(final ConsensusReads consensusReads)
    {
        try
        {
            ConsensusReadInfo consensusReadInfo = consensusReads.createConsensusRead(mReads, mFragmentCoords, mUmi);

            // set consensus read attributes
            ConsensusType consensusType = mDualStrand ? DUAL : SINGLE;

            int pcrClusterCount = isIllumina() ? calculatePCRClusterCount(this) : UNSET_COUNT;
            int firstInPairCount = calculateFirstInPairCount();

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, totalReadCount(), firstInPairCount, consensusType, pcrClusterCount);

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

    private int calculateFirstInPairCount()
    {
        // poly-G reads do not count towards the first in pair count, but collapsed UMI group reads do
        int firstInPairCount = (int)mReads.stream().filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag()).count();
        firstInPairCount += (int)mNonConsensusReads.stream().filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag()).count();

        int totalReadCount = mReads.size() + mNonConsensusReads.size();
        boolean isPrimaryGroup = firstInPairCount >= 0.5f * totalReadCount;

        if(!isPrimaryGroup)
            firstInPairCount = totalReadCount - firstInPairCount; // adjusted so both reads report the same ratio

        return firstInPairCount;
    }

    public String toString()
    {
        return format("id(%s) reads(%d) coords(%s)", mUmi, totalReadCount(), mFragmentCoords.Key);
    }
}
