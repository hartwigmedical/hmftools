package com.hartwig.hmftools.patientreporter.slicing;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class Slicer {

    @NotNull
    private final List<GenomeRegion> regions;

    public Slicer(@NotNull final List<GenomeRegion> regions) {
        this.regions = regions;
    }

    @VisibleForTesting
    int numberOfRegions() {
        return regions.size();
    }

    @VisibleForTesting
    long numberOfBases() {
        long bases = 0;
        for (GenomeRegion region : regions) {
            bases += region.bases();
        }
        return bases;
    }
}
