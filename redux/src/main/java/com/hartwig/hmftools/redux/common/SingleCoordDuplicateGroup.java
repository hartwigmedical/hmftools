package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.List;
import java.util.function.BiConsumer;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class SingleCoordDuplicateGroup extends DuplicateGroup
{
    private final List<SAMRecord> mReads;
    private final FragmentCoords mFragmentCoords;

    public SingleCoordDuplicateGroup(final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        this(null, reads, fragmentCoords);
    }

    public SingleCoordDuplicateGroup(final String id, final List<SAMRecord> reads, final FragmentCoords fragmentCoords)
    {
        super(id);

        mReads = reads;
        mFragmentCoords = fragmentCoords;

    }

    public void addRead(final SAMRecord read)
    {
        mReads.add(read);
    }

    public void addReads(final List<SAMRecord> reads)
    {
        mReads.addAll(reads);
    }

    @Override
    public List<SAMRecord> reads()
    {
        return mReads;
    }

    @Override
    public boolean readIsLower()
    {
        return mFragmentCoords.ReadIsLower;
    }

    @Override
    public int readCount()
    {
        return mReads.size();
    }

    public FragmentCoords fragmentCoordinates()
    {
        return mFragmentCoords;
    }

    public List<SAMRecord> duplicate()
    {
        return mReads;
    }

    public String coordinatesKey()
    {
        return mFragmentCoords.Key;
    }

    @Override
    public String consensusCoordinatesKey()
    {
        return mFragmentCoords.Key;
    }

    @Override
    public void processReads(final BiConsumer<SAMRecord, FragmentCoords> consumer)
    {
        mReads.forEach(read -> consumer.accept(read, mFragmentCoords));
    }

    public String toString()
    {
        return format("id(%s) reads(%d) coords(%s)", mUmiId, mReads.size(), mFragmentCoords.Key);
    }
}
