package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LocalPhaseSetTest {

    @Test
    public void testRightInLeftDel() {
        VariantHotspot left = create(100, "ATT", "T");
        assertTrue(LocalPhaseSet.rightInLeftDel(left, create(101, "T", "C")));
        assertTrue(LocalPhaseSet.rightInLeftDel(left, create(102, "T", "C")));
        assertFalse(LocalPhaseSet.rightInLeftDel(left, create(103, "T", "C")));
    }

    @Test
    public void testOffsetSNVs() {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "C", "T");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));
    }

    @Test
    public void testDelToTheRight() {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "CAT", "T");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));
    }

    @Test
    public void testInsToTheRight() {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "C", "TAT");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));
    }

    @Test
    public void testDelToTheLeft() {
        VariantHotspot left = create(100, "AA", "A");
        VariantHotspot right = create(102, "G", "T");
        assertEquals(-1, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "AA", "A");
        right = create(103, "G", "T");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "AAC", "A");
        right = create(103, "G", "T");
        assertEquals(-1, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "AAC", "A");
        right = create(104, "G", "T");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));
    }

    @Test
    public void testInsToTheLeft() {
        VariantHotspot left = create(100, "A", "AA");
        VariantHotspot right = create(101, "G", "T");
        assertEquals(-2, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "A", "AA");
        right = create(102, "G", "T");
        assertEquals(-3, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "A", "AAC");
        right = create(101, "G", "T");
        assertEquals(-3, LocalPhaseSet.adjustedOffset(left, right));

        left = create(100, "A", "AAC");
        right = create(102, "G", "T");
        assertEquals(-4, LocalPhaseSet.adjustedOffset(left, right));
    }

    @Test
    public void testTwoInsertsAtSameLocation() {
        VariantHotspot left = create(100, "A", "AAA");
        VariantHotspot right = create(100, "A", "AA");
        assertEquals(0, LocalPhaseSet.adjustedOffset(left, right));
        assertEquals(0, LocalPhaseSet.adjustedOffset(right, left));
    }

    @Test
    public void testTwoDeletesAtSameLocation() {
        VariantHotspot left = create(100, "AAA", "A");
        VariantHotspot right = create(100, "AA", "A");
        assertEquals(0, LocalPhaseSet.adjustedOffset(left, right));
        assertEquals(0, LocalPhaseSet.adjustedOffset(right, left));
    }

    @NotNull
    static VariantHotspot create(long position, String ref, String alt) {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").ref(ref).alt(alt).position(position).build();
    }

}
