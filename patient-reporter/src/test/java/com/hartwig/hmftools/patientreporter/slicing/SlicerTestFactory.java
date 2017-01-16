package com.hartwig.hmftools.patientreporter.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;

import org.jetbrains.annotations.NotNull;

public final class SlicerTestFactory {

    private SlicerTestFactory() {
    }

    @NotNull
    public static Slicer oneRegionTestSlicer(@NotNull final String chromosome, final long start, final long end) {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(chromosome, new GenomeRegion("1", start, end));
        return new Slicer(regionMap);
    }
}
