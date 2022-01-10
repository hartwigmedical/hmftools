package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.common.EvictingArray.calculateSize;
import static com.hartwig.hmftools.sage.common.EvictingArray.calculateSizeOld;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.candidate.RefContext;

import org.junit.Before;
import org.junit.Test;

public class EvictingArrayTest
{
    private static final int CAPACITY = 256;

    private EvictingArray mArray;
    private EvictionHandler mHandler;

    @Before
    public void setup()
    {
        mHandler = new EvictionHandler();
        mArray = new EvictingArray(CAPACITY, mHandler);
    }

    @Test
    public void testCapacityCalc()
    {
        assertEquals(8, calculateSizeOld(8));
        assertEquals(8, calculateSize(8));

        assertEquals(32, calculateSizeOld(17));
        assertEquals(32, calculateSize(32));

        assertEquals(256, calculateSizeOld(129));
        assertEquals(256, calculateSize(256));

        assertEquals(1024, calculateSizeOld(513));
        assertEquals(1024, calculateSize(1024));

        assertEquals(8192, calculateSizeOld(4097));
        assertEquals(8192, calculateSize(8192));
    }

    @Test
    public void testInitialPosition()
    {
        mArray.computeIfAbsent(512, EvictingArrayTest::create);
        assertEquals(257, mArray.minPosition());
    }

    @Test
    public void testFillCapacity()
    {
        for(int i = 0; i < CAPACITY; i++)
        {
            mArray.computeIfAbsent(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(1000, mArray.minPosition());
        assertEquals(0, mHandler.list.size());
    }

    @Test
    public void testSelectFromExistingDoesNotEvict()
    {
        testFillCapacity();
        testFillCapacity();
    }

    @Test
    public void testCapacityOverflow()
    {
        for(int i = 0; i < CAPACITY + 100; i++)
        {
            mArray.computeIfAbsent(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(1100, mArray.minPosition());
        assertEquals(100, mHandler.list.size());
        assertEquals(mHandler.list.get(0).position(), 1000);
        assertEquals(mHandler.list.get(99).position(), 1099);
    }

    static class EvictionHandler implements Consumer<RefContext>
    {
        private final List<GenomePosition> list = Lists.newArrayList();

        @Override
        public void accept(final RefContext position)
        {
            list.add(position);
        }
    }

    private static RefContext create(int pos)
    {
        return new RefContext("1", pos, false);
    }
}
