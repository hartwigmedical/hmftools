package com.hartwig.hmftools.patientreporter.slicing;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class Slicer {

    @NotNull
    private final Multimap<String, GenomeRegion> regions;

    Slicer(@NotNull final Multimap<String, GenomeRegion> regions) {
        this.regions = regions;
    }

    @VisibleForTesting
    int numberOfRegions() {
        return regions.size();
    }

    @VisibleForTesting
    long numberOfBases() {
        long bases = 0;
        for (GenomeRegion region : regions.values()) {
            bases += region.bases();
        }
        return bases;
    }
}
