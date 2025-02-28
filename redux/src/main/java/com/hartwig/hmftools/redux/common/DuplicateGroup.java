package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNSET_COUNT;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.UmiReadType.DUAL;
import static com.hartwig.hmftools.common.bam.UmiReadType.SINGLE;
import static com.hartwig.hmftools.common.collect.Cluster.clusterCount;
import static com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.getReadNameAttributes;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.IlluminaReadNameAttributes;
import com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.TileCoord;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final String mUmiId; // the UMI if enabled

    // with duplicate group collapsing some reads in mReads may not have mFragmentCoords FragmentCoords
    private final FragmentCoords mFragmentCoords;
    private final List<SAMRecord> mReads;

    private SAMRecord mConsensusRead;
    private SAMRecord mPrimaryRead; // if no consensus is formed, the selected primary read
    private boolean mDualStrand;

    private int mPCRClusterCount;

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

        mPCRClusterCount = UNSET_COUNT;
    }

    public void addRead(final SAMRecord read) { mReads.add(read); }
    public void addReads(final List<SAMRecord> reads) { mReads.addAll(reads); }

    public List<SAMRecord> reads() { return mReads; }
    public int readCount() { return mReads.size(); }

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
            int firstInPairCount = (int)mReads.stream().filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag()).count();
            int readCount = mReads.size();
            boolean isPrimaryGroup = firstInPairCount >= readCount / 2;

            if(!isPrimaryGroup)
                firstInPairCount = readCount - firstInPairCount; // adjusted so both reads report the same ratio

            UmiReadType umiReadType = mDualStrand ? DUAL : SINGLE;

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount, firstInPairCount, umiReadType, mPCRClusterCount);

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

    private static final Set<Integer> OPTICAL_DUPLICATE_TILE_DIFFERENCE = Sets.newHashSet(0, 1, 999, 1000, 1001);
    private static final int OPTICAL_DUPLICATE_DISTANCE_THRESHOLD = 2_500;

    public void setPCRClusterCount(final SequencingType sequencingType)
    {
        if(sequencingType != ILLUMINA)
            return;

        List<IlluminaReadNameAttributes> readNameAttributes = Lists.newArrayList();
        for(SAMRecord read : mReads)
        {
            IlluminaReadNameAttributes attributes = getReadNameAttributes(read.getReadName());
            if(attributes == null)
                return;

            readNameAttributes.add(attributes);
        }

        if(mReads.size() == 2)
        {
            IlluminaReadNameAttributes readNameAttributes1 = readNameAttributes.get(0);
            IlluminaReadNameAttributes readNameAttributes2 = readNameAttributes.get(1);
            int tileDifference = abs(readNameAttributes1.tileNumber() - readNameAttributes2.tileNumber());
            if(readNameAttributes1.laneKey().equals(readNameAttributes2.laneKey())
                    && OPTICAL_DUPLICATE_TILE_DIFFERENCE.contains(tileDifference))
            {
                mPCRClusterCount = 1;
                return;
            }

            mPCRClusterCount = 2;
            return;
        }

        Map<String, List<TileCoord>> tileCoordsByTile = Maps.newHashMap();
        for(IlluminaReadNameAttributes attributes : readNameAttributes)
        {
            String tileKey = attributes.tileKey();
            TileCoord tileCoord = attributes.tileCoord();
            tileCoordsByTile.computeIfAbsent(tileKey, key -> Lists.newArrayList());
            tileCoordsByTile.get(tileKey).add(tileCoord);
        }

        mPCRClusterCount = 0;
        for(List<TileCoord> tileCoords : tileCoordsByTile.values())
        {
            if(tileCoords.size() == 1)
            {
                mPCRClusterCount++;
                continue;
            }

            if(tileCoords.size() == 2)
            {
                TileCoord tileCoord1 = tileCoords.get(0);
                TileCoord tileCoord2 = tileCoords.get(1);
                if(tileCoord1.distance(tileCoord2) <= OPTICAL_DUPLICATE_DISTANCE_THRESHOLD)
                    mPCRClusterCount++;
                else
                    mPCRClusterCount += 2;

                continue;
            }

            mPCRClusterCount += clusterCount(tileCoords, TileCoord::distance, OPTICAL_DUPLICATE_DISTANCE_THRESHOLD);
        }
    }

    public String toString()
    {
        return format("id(%s) reads(%d) coords(%s)", mUmiId, mReads.size(), mFragmentCoords.Key);
    }
}
