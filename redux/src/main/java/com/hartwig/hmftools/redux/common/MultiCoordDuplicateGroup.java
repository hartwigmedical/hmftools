package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class MultiCoordDuplicateGroup extends DuplicateGroup
{
    private final List<ReadWithFragCoords> mReads;

    public MultiCoordDuplicateGroup(final List<ReadWithFragCoords> reads)
    {
        this(null, reads);
    }

    public MultiCoordDuplicateGroup(final String id, final List<ReadWithFragCoords> reads)
    {
        super(id);

        mReads = reads;
    }

    @Override
    public int readCount()
    {
        return mReads.size();
    }

    @Override
    public List<SAMRecord> reads()
    {
        return mReads.stream().map(x -> x.Read).collect(Collectors.toList());
    }

    public Set<String> coordinateKeys()
    {
        return mReads.stream().map(x -> x.FragCoords.Key).collect(Collectors.toCollection(() -> Sets.newTreeSet()));
    }

    @Override
    public boolean readIsLower()
    {
        return mReads.get(0).FragCoords.ReadIsLower;
    }

    @Override
    public String consensusCoordinatesKey()
    {
        return FragmentCoords.fromRead(mConsensusRead, false).Key;
    }

    @Override
    public void processReads(final BiConsumer<SAMRecord, FragmentCoords> consumer)
    {
        mReads.forEach(x -> consumer.accept(x.Read, x.FragCoords));
    }

    public String toString()
    {
        String coordinateKeysStr = coordinateKeys().stream().collect(Collectors.joining(";"));
        return format("id(%s) reads(%d) coords(%s)", mUmiId, mReads.size(), coordinateKeysStr);
    }
}
