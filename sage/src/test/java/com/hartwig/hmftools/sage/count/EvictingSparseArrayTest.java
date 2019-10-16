package com.hartwig.hmftools.sage.count;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.junit.Before;
import org.junit.Test;

public class EvictingSparseArrayTest {

    private EvictingSparseArray<GenomePosition> victim;
    private EvictionHandler hander;

    @Before
    public void setup() {
        hander = new EvictionHandler();
        victim = new EvictingSparseArray<>(256, hander);
    }

    @Test
    public void testStuff() {

        victim.computeIfAbsent(1001, EvictingSparseArrayTest::create);
        assertEquals(746, victim.minPosition());

        victim.computeIfAbsent(1001, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1002, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1005, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1100, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1256, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1300, EvictingSparseArrayTest::create);
        victim.computeIfAbsent(1256, EvictingSparseArrayTest::create);


        System.out.println("sdf");
    }


    class EvictionHandler implements Consumer<GenomePosition> {

        private final List<GenomePosition> list = Lists.newArrayList();

        @Override
        public void accept(final GenomePosition position) {
            list.add(position);
        }
    }

    static GenomePosition create(long pos) {
        return GenomePositions.create("CHROM", pos);
    }

}
