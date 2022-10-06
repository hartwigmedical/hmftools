package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.common.EvictingArray.MIN_CAPACITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.candidate.RefContext;

import org.junit.Before;
import org.junit.Test;

public class EvictingArrayTest
{
    private EvictingArray mArray;
    private EvictionHandler mHandler;

    private static final int TEST_CAPACITY = 100;

    @Before
    public void setup()
    {
        mHandler = new EvictionHandler();
        mArray = new EvictingArray(TEST_CAPACITY, mHandler);
    }

    @Test
    public void testDepth()
    {
        // no initial depth
        assertNull(mArray.getDepth(100));
        assertNull(mArray.getDepth(200));

        // first check depth only
        mArray.registerDepth(100);
        assertEquals(50, mArray.minPosition());
        mArray.registerDepth(50);
        mArray.registerDepth(75);
        mArray.registerDepth(100);
        mArray.registerDepth(149);
        assertEquals(50, mArray.minPosition());
        assertEquals(1, mArray.getDepth(50).intValue());
        assertEquals(1, mArray.getDepth(75).intValue());
        assertEquals(2, mArray.getDepth(100).intValue());
        assertEquals(0, mArray.getDepth(148).intValue()); // checking boundaries
        assertEquals(1, mArray.getDepth(149).intValue());
        assertNull(mArray.getDepth(150)); // past the end

        // shift so that the earliest position has been flushed
        mArray.registerDepth(150);
        assertEquals(51, mArray.minPosition());
        assertNull(mArray.getDepth(50));
        assertEquals(1, mArray.getDepth(75).intValue());
        assertEquals(2, mArray.getDepth(100).intValue());
        assertEquals(1, mArray.getDepth(149).intValue());
        assertEquals(1, mArray.getDepth(150).intValue());

        // shift again
        mArray.registerDepth(180);
        assertEquals(81, mArray.minPosition());
        assertNull(mArray.getDepth(75));
        assertEquals(2, mArray.getDepth(100).intValue());
        assertEquals(1, mArray.getDepth(149).intValue());
        assertEquals(1, mArray.getDepth(180).intValue());

        // correct indices are maintained
        mArray.registerDepth(100);
        mArray.registerDepth(180);

        assertEquals(3, mArray.getDepth(100).intValue());
        assertEquals(2, mArray.getDepth(180).intValue());

        // move completely beyond the current range
        mArray.registerDepth(300);
        mArray.registerDepth(350);

        assertNull(mArray.getDepth(180));

        for(int pos = mArray.minPosition(); pos < mArray.minPosition() + mArray.capacity(); ++pos)
        {
            int expectDepth = pos == 300 || pos == 350 ? 1 : 0;
            assertEquals(expectDepth, mArray.getDepth(pos).intValue());
        }
    }

    @Test
    public void testFillCapacity()
    {
        for(int i = 0; i < TEST_CAPACITY; i++)
        {
            mArray.getOrCreateRefContext(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(0, mHandler.items().size());
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
        for(int i = 0; i < TEST_CAPACITY + 100; i++)
        {
            mArray.getOrCreateRefContext(1000 + i, EvictingArrayTest::create);
        }

        assertEquals(1100, mArray.minPosition());
        assertEquals(100, mHandler.items().size());
        assertEquals(1000, mHandler.items().get(0).position());
        assertEquals(1099, mHandler.items().get(99).position());
    }

    static class EvictionHandler implements Consumer<RefContext>
    {
        private final List<GenomePosition> mItems = Lists.newArrayList();

        public List<GenomePosition> items() { return mItems; }

        @Override
        public void accept(final RefContext position)
        {
            mItems.add(position);
        }
    }

    private static RefContext create(int pos)
    {
        return new RefContext("1", pos);
    }
}
