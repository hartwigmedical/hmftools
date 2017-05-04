package com.hartwig.hmftools.common.slicing;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public enum SlicerFactory {
    ;

    @NotNull
    public static Slicer fromBedFile(@NotNull String bedFile) throws IOException, EmptyFileException {
        return new UnsortedSlicer(BEDFileLoader.fromBedFile(bedFile));
    }

    @NotNull
    public static Slicer sortedSlicer(@NotNull String bedFile) throws IOException, EmptyFileException {
        return new SortedSlicer(BEDFileLoader.fromBedFile(bedFile));
    }

    @NotNull
    public static Slicer fromSingleGenomeRegion(@NotNull final GenomeRegion region) {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(region.chromosome(), region);
        return new UnsortedSlicer(regionMap);
    }

}
