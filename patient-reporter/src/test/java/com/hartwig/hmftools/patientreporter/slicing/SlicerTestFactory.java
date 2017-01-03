package com.hartwig.hmftools.patientreporter.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;

import org.jetbrains.annotations.NotNull;

public final class SlicerTestFactory {

    private SlicerTestFactory() {
    }

    @NotNull
    public static Slicer oneRegionTestSlicer(@NotNull String chromosome, long start, long end) {
        SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(chromosome, new GenomeRegion(start, end));
        return new Slicer(regionMap);
    }
}
