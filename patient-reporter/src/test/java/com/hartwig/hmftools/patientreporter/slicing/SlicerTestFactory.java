package com.hartwig.hmftools.patientreporter.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;

import org.jetbrains.annotations.NotNull;

public final class SlicerTestFactory {

    private SlicerTestFactory() {
    }

    @NotNull
    public static Slicer forGenomeRegion(@NotNull final GenomeRegion region) {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(region.chromosome(), region);
        return new Slicer(regionMap);
    }
}
