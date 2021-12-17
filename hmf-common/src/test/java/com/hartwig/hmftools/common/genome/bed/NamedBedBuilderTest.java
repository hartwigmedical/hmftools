package com.hartwig.hmftools.common.genome.bed;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class NamedBedBuilderTest {

    @Test
    public void testNoOverlap() {
        final NamedBedBuilder victim = new NamedBedBuilder();

        final List<NamedBed> orderedInput = Lists.newArrayList();
        orderedInput.add(create(101, 200, "A"));
        orderedInput.add(create(301, 400, "B"));
        orderedInput.add(create(501, 600, "C"));

        final List<NamedBed> randomInput = new ArrayList<>(orderedInput);
        Collections.shuffle(randomInput);
        randomInput.forEach(victim::addBed);

        final List<NamedBed> result = victim.build();
        for (int i = 0; i < orderedInput.size(); i++) {
            assertEquals(orderedInput.get(i), result.get(i));
        }
    }

    @Test
    public void testFullyContains() {
        final NamedBedBuilder victim = new NamedBedBuilder();

        final List<NamedBed> orderedInput = Lists.newArrayList();
        orderedInput.add(create(101, 200, "A"));
        orderedInput.add(create(301, 400, "B"));
        orderedInput.add(create(501, 600, "C"));
        orderedInput.forEach(victim::addBed);

        assertFalse(victim.addBed(create(150, 190, "A")));
        assertFalse(victim.addBed(create(150, 190, "SXD")));
        assertFalse(victim.addBed(create(301, 400, "SXD")));
        assertFalse(victim.addBed(create(599, 600, "SXD")));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testPartialOverlap() {
        final NamedBedBuilder victim = new NamedBedBuilder();

        final List<NamedBed> orderedInput = Lists.newArrayList();
        orderedInput.add(create(101, 200, "A"));
        orderedInput.add(create(301, 400, "B"));
        orderedInput.add(create(501, 600, "C"));
        orderedInput.forEach(victim::addBed);

        assertFalse(victim.addBed(create(40, 101, "A")));
    }

    @NotNull
    private static NamedBed create(final int start, final int end, final String name) {
        return ImmutableNamedBed.builder().chromosome("1").start(start).end(end).name(name).build();
    }
}
