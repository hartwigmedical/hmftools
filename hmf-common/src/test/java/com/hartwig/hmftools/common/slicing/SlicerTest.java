package com.hartwig.hmftools.common.slicing;

import static org.junit.Assert.assertSame;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantImpl;

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

    private static void excludedChromosomes(@NotNull Slicer slicer) {
        assertSlicerInclude(false, buildVariant("1", 0), slicer);
        assertSlicerInclude(false, buildVariant("NotExists", 100), slicer);
    }

    private static void sortedVariants(@NotNull Slicer slicer) {
        assertSlicerInclude(false, buildVariant("X", 1), slicer);
        assertSlicerInclude(true, buildVariant("X", 150), slicer);
        assertSlicerInclude(false, buildVariant("X", 250), slicer);
        assertSlicerInclude(true, buildVariant("X", 300), slicer);
        assertSlicerInclude(true, buildVariant("X", 400), slicer);
        assertSlicerInclude(true, buildVariant("Y", 570), slicer);
        assertSlicerInclude(false, buildVariant("Y", 1000), slicer);
        assertSlicerInclude(false, buildVariant("NotExists", 100), slicer);
    }

    private static void unsortedVariants(@NotNull Slicer slicer) {
        assertSlicerInclude(true, buildVariant("X", 150), slicer);
        assertSlicerInclude(true, buildVariant("X", 300), slicer);
        assertSlicerInclude(true, buildVariant("X", 400), slicer);
        assertSlicerInclude(true, buildVariant("Y", 570), slicer);
        assertSlicerInclude(false, buildVariant("X", 1), slicer);
        assertSlicerInclude(false, buildVariant("X", 250), slicer);
        assertSlicerInclude(false, buildVariant("Y", 1000), slicer);
        assertSlicerInclude(false, buildVariant("NotExists", 100), slicer);
    }

    @NotNull
    private static SomaticVariant buildVariant(@NotNull String chromosome, long position) {
        return new SomaticVariantImpl.Builder().chromosome(chromosome).position(position).build();
    }

    private static void assertSlicerInclude(boolean expectedResult, @NotNull SomaticVariant variant, @NotNull Slicer slicer) {
        assertSame(expectedResult, slicer.includes(variant));
    }
}