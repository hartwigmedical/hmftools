package com.hartwig.hmftools.common.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import com.hartwig.hmftools.common.variant.Variant;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SlicerTest {

    private Slicer unsortedSlicer;
    private Slicer sortedSlicer;

    @Before
    public void setup() {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put("X", new GenomeRegion("X", 100, 200));
        regionMap.put("X", new GenomeRegion("X", 300, 400));
        regionMap.put("Y", new GenomeRegion("Y", 500, 600));

        unsortedSlicer = new UnsortedSlicer(regionMap);
        sortedSlicer = new SortedSlicer(regionMap);
    }

    @Test
    public void excludedChromosomes() {
        excludedChromosomes(sortedSlicer);
        excludedChromosomes(unsortedSlicer);
    }

    @Test
    public void sortedVariants() {
        sortedVariants(unsortedSlicer);
        sortedVariants(sortedSlicer);
    }

    @Test
    public void unsortedSlicerWithUnsortedVariants() {
        unsortedVariants(unsortedSlicer);
    }

    @Test(expected = IllegalArgumentException.class)
    public void sortedSlicerWithUnsortedVariants() {
        unsortedVariants(sortedSlicer);
    }

    private void excludedChromosomes(Slicer aSlicer) {
        assertSlicers(false, buildVariant("1", 0), aSlicer);
        assertSlicers(false, buildVariant("NotExists", 100), aSlicer);
    }

    private void sortedVariants(Slicer aSlicer) {
        assertSlicers(false, buildVariant("X", 1), aSlicer);
        assertSlicers(true, buildVariant("X", 150), aSlicer);
        assertSlicers(false, buildVariant("X", 250), aSlicer);
        assertSlicers(true, buildVariant("X", 300), aSlicer);
        assertSlicers(true, buildVariant("X", 400), aSlicer);
        assertSlicers(true, buildVariant("Y", 570), aSlicer);
        assertSlicers(false, buildVariant("Y", 1000), aSlicer);
        assertSlicers(false, buildVariant("NotExists", 100), aSlicer);
    }

    private void unsortedVariants(Slicer aSlicer) {
        assertSlicers(true, buildVariant("X", 150), aSlicer);
        assertSlicers(true, buildVariant("X", 300), aSlicer);
        assertSlicers(true, buildVariant("X", 400), aSlicer);
        assertSlicers(true, buildVariant("Y", 570), aSlicer);
        assertSlicers(false, buildVariant("X", 1), aSlicer);
        assertSlicers(false, buildVariant("X", 250), aSlicer);
        assertSlicers(false, buildVariant("Y", 1000), aSlicer);
        assertSlicers(false, buildVariant("NotExists", 100), aSlicer);
    }

    @NotNull
    private static SomaticVariant buildVariant(@NotNull String chromosome, long position) {
        return new SomaticVariant.Builder().chromosome(chromosome).position(position).build();
    }

    private void assertSlicers(boolean expectedResult, Variant variant, Slicer aSlicer) {
        assertSame(expectedResult, aSlicer.includes(variant));
    }
}