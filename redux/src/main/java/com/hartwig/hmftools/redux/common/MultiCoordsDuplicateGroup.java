package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.UmiReadType.SINGLE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class MultiCoordsDuplicateGroup
{
    // NOTE: We assume that all the FragmentCoords of added reads have the same ReadIsLower value.

    private final Map<FragmentCoords, List<SAMRecord>> mReadsByCoords;

    private SAMRecord mConsensusRead;
    private SAMRecord mPrimaryRead; // if no consensus is formed, the selected primary read

    public MultiCoordsDuplicateGroup(final Map<FragmentCoords, List<SAMRecord>> readsByCoords)
    {
        mReadsByCoords = readsByCoords;

        mConsensusRead = null;
        mPrimaryRead = null;
    }

    public MultiCoordsDuplicateGroup(final List<SAMRecord> reads, final FragmentCoords fragCoords)
    {
        this(Maps.newHashMap());

        mReadsByCoords.put(fragCoords, Lists.newArrayList(reads));
    }

    public MultiCoordsDuplicateGroup(final SAMRecord read, final FragmentCoords coords)
    {
        this(Lists.newArrayList(read), coords);
    }

    public MultiCoordsDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        this(duplicateGroup.reads(), duplicateGroup.fragmentCoordinates());
    }

    public void addRead(final SAMRecord read, final FragmentCoords fragmentCoords)
    {
        List<SAMRecord> reads = mReadsByCoords.get(fragmentCoords);
        if(reads == null)
        {
            reads = Lists.newArrayList();
            mReadsByCoords.put(fragmentCoords, reads);
        }

        reads.add(read);
    }

    public void addReads(final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        List<SAMRecord> currentReads = mReadsByCoords.get(fragmentCoords);
        if(currentReads == null)
        {
            currentReads = Lists.newArrayList();
            mReadsByCoords.put(fragmentCoords, currentReads);
        }

        currentReads.addAll(reads);
    }

    public MultiCoordsDuplicateGroup merge(final MultiCoordsDuplicateGroup group)
    {
        for(Map.Entry<FragmentCoords, List<SAMRecord>> coordsAndReads : group.mReadsByCoords.entrySet())
        {
            FragmentCoords fragCoords = coordsAndReads.getKey();
            List<SAMRecord> reads = coordsAndReads.getValue();
            addReads(reads, fragCoords);
        }

        return this;
    }

    public MultiCoordsDuplicateGroup merge(final DuplicateGroup group)
    {
        addReads(group.reads(), group.fragmentCoordinates());
        return this;
    }

    public int readCount()
    {
        return mReadsByCoords.values().stream().mapToInt(List::size).sum();
    }

    public int coordCount()
    {
        return mReadsByCoords.size();
    }

    public Map<FragmentCoords, List<SAMRecord>> readsByCoords()
    {
        return mReadsByCoords;
    }

    public List<SAMRecord> reads()
    {
        if(mReadsByCoords.size() == 1)
            return mReadsByCoords.values().stream().findAny().orElse(null);

        List<SAMRecord> reads = Lists.newArrayList();
        mReadsByCoords.values().forEach(reads::addAll);
        return reads;
    }

    public boolean readIsLower()
    {
        return anyFragmentCoordinates().ReadIsLower;
    }

    public SAMRecord consensusRead()
    {
        return mConsensusRead;
    }

    public void setPrimaryRead(final SAMRecord read)
    {
        mPrimaryRead = read;
    }

    public boolean isPrimaryRead(final SAMRecord read)
    {
        return mPrimaryRead == read;
    }

    public FragmentCoords anyFragmentCoordinates()
    {
        return mReadsByCoords.keySet().stream().findAny().orElse(null);
    }

    public String coordinateKeysStr()
    {
        return mReadsByCoords.keySet().stream().map(x -> x.Key).sorted().collect(Collectors.joining(";"));
    }

    public String consensusCoordinatesKey()
    {
        return FragmentCoords.fromRead(mConsensusRead, false).Key;
    }

    public boolean isSingleRead()
    {
        if(mReadsByCoords.size() > 1)
            return false;

        return readCount() == 1;
    }

    public boolean isSingleCoordDuplicateGroup()
    {
        return mReadsByCoords.size() == 1;
    }

    public ReadInfo toReadInfo()
    {
        if(!isSingleRead())
            throw new IllegalStateException("Cannot create a ReadInfo from a DuplicateGroup with more than one read");

        for(Map.Entry<FragmentCoords, List<SAMRecord>> coordsAndReads : mReadsByCoords.entrySet())
            return new ReadInfo(coordsAndReads.getValue().get(0), coordsAndReads.getKey());

        throw new RuntimeException("Unreachable");
    }

    public DuplicateGroup toSingleCoordDuplicateGroup()
    {
        if(!isSingleCoordDuplicateGroup())
            throw new IllegalStateException("Cannot create a DuplicateGroup from a MultiCoordDuplicateGroup without unique FragmentCoords");

        for(Map.Entry<FragmentCoords, List<SAMRecord>> coordsAndReads : mReadsByCoords.entrySet())
            return new DuplicateGroup(coordsAndReads.getValue(), coordsAndReads.getKey());

        throw new RuntimeException("Unreachable");
    }

    public void formConsensusRead(final ConsensusReads consensusReads)
    {
        List<SAMRecord> reads = reads();
        boolean readIsLower = readIsLower();

        try
        {
            ConsensusReadInfo consensusReadInfo = consensusReads.createConsensusRead(reads, readIsLower, null);

            // set consensus read attributes
            int firstInPairCount = (int) reads.stream().filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag()).count();
            int readCount = reads.size();
            boolean isPrimaryGroup = firstInPairCount >= readCount / 2;

            if(!isPrimaryGroup)
                firstInPairCount = readCount - firstInPairCount; // adjusted so both reads report the same ratio

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount, firstInPairCount, SINGLE);

            mConsensusRead = consensusReadInfo.ConsensusRead;
        }
        catch(Exception e)
        {
            RD_LOGGER.error("error forming consensus: {}", toString());

            for(SAMRecord read : reads)
                RD_LOGGER.error("read: {}", readToString(read));

            e.printStackTrace();
        }
    }

    @Override
    public String toString()
    {
        return format("reads(%d) coords(%s)", readCount(), coordinateKeysStr());
    }
}
