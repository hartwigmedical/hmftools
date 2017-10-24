package com.hartwig.hmftools.common.slicing;

import static org.junit.Assert.assertSame;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class SlicerTest {

    private Slicer forwardSlicer;
    private Slicer bidirectionalSlicer;

    @Before
    public void setup() {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put("X", GenomeRegionFactory.create("X", 100, 200));
        regionMap.put("X", GenomeRegionFactory.create("X", 300, 400));
        regionMap.put("Y", GenomeRegionFactory.create("Y", 500, 600));

        bidirectionalSlicer = new BidirectionalSlicer(regionMap);
        forwardSlicer = new ForwardSlicer(regionMap);
    }

    @Test
    public void excludedChromosomes() {
        excludedChromosomes(forwardSlicer);
        excludedChromosomes(bidirectionalSlicer);
    }

    @Test
    public void sortedVariants() {
        sortedVariants(bidirectionalSlicer);
        sortedVariants(forwardSlicer);
    }

    @Test
    public void bidirectionalSlicerWithUnsortedVariants() {
        unsortedVariants(bidirectionalSlicer);
    }

    @Test(expected = IllegalArgumentException.class)
    public void forwardSlicerWithUnsortedVariants() {
        unsortedVariants(forwardSlicer);
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