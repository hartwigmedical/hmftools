package com.hartwig.hmftools.common.slicing;

import java.io.IOException;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;

import org.jetbrains.annotations.NotNull;

public enum SlicerFactory {
    ;

    @NotNull
    public static Slicer fromBedFile(@NotNull String bedFile) throws IOException {
        return new BidirectionalSlicer(BEDFileLoader.fromBedFile(bedFile));
    }

    @NotNull
    public static Slicer fromRegions(@NotNull final SortedSetMultimap<String, ? extends GenomeRegion> regions) {
        return new BidirectionalSlicer(regions);
    }
}
