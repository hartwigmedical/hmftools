package com.hartwig.hmftools.sage.count;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.junit.Before;
import org.junit.Test;

public class EvictingArrayTest
{
    private static final int CAPACITY = 256;

    private EvictingArray<GenomePosition> victim;
    private EvictionHandler handler;

    @Before
    public void setup()
    {
        handler = new EvictionHandler();
        victim = new EvictingArray<>(CAPACITY, handler);
    }

    @Test
    public void testCapacity()
    {
        victim = new EvictingArray<>(151, handler);
        assertEquals(256, victim.capacity());

        victim = new EvictingArray<>(CAPACITY, handler);
        assertEquals(CAPACITY, victim.capacity());

        victim = new EvictingArray<>(300, handler);
        assertEquals(512, victim.capacity());
    }

    @Test
    public void testInitialPosition()
    {
        victim.computeIfAbsent(512, EvictingArrayTest::create);
        assertEquals(257, victim.minPosition());
    }

    @Test
    public void testFillCapacity()
    {
        for(int i = 0; i < CAPACITY; i++)
        {
            victim.computeIfAbsent(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(1000, victim.minPosition());
        assertEquals(0, handler.list.size());
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
            victim.computeIfAbsent(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(1100, victim.minPosition());
        assertEquals(100, handler.list.size());
        assertEquals(handler.list.get(0).position(), 1000);
        assertEquals(handler.list.get(99).position(), 1099);
    }

    static class EvictionHandler implements Consumer<GenomePosition>
    {
        private final List<GenomePosition> list = Lists.newArrayList();

        @Override
        public void accept(final GenomePosition position)
        {
            list.add(position);
        }
    }

    private static GenomePosition create(long pos)
    {
        return GenomePositions.create("CHROM", pos);
    }
}
