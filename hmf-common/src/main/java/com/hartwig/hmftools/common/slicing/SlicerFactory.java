package com.hartwig.hmftools.common.slicing;

import java.io.IOException;

import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfSlicerFileLoader;

import org.jetbrains.annotations.NotNull;

public enum SlicerFactory {
    ;

    @NotNull
    public static Slicer fromBedFile(@NotNull String bedFile) throws IOException, EmptyFileException {
        return new BidirectionalSlicer(BEDFileLoader.fromBedFile(bedFile));
    }

    @NotNull
    public static HmfSlicer fromHmfSlicerFile(@NotNull String hmfSlicer) throws IOException, EmptyFileException {
        return new HmfSlicer(HmfSlicerFileLoader.fromHmfSlicerFile(hmfSlicer));
    }

    @NotNull
    public static Slicer forwardSlicer(@NotNull String bedFile) throws IOException, EmptyFileException {
        return new ForwardSlicer(BEDFileLoader.fromBedFile(bedFile));
    }

    @NotNull
    public static Slicer fromSingleGenomeRegion(@NotNull final GenomeRegion region) {
        final SortedSetMultimap<String, GenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(region.chromosome(), region);
        return new BidirectionalSlicer(regionMap);
    }

    @NotNull
    public static HmfSlicer hmfSlicerFromSingleGenomeRegion(@NotNull final HmfGenomeRegion region) {
        final SortedSetMultimap<String, HmfGenomeRegion> regionMap = TreeMultimap.create();
        regionMap.put(region.chromosome(), region);
        return new HmfSlicer(regionMap);
    }
}
