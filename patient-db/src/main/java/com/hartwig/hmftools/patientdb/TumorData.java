package com.hartwig.hmftools.patientdb;

import java.util.List;

import org.jetbrains.annotations.Nullable;

class TumorData {
    @Nullable
    private final String location;
    @Nullable
    private final List<String> biopsyLocations;
    @Nullable
    private final String entryStage;

    TumorData(@Nullable final String location, @Nullable final List<String> biopsyLocations,
            @Nullable final String entryStage) {
        this.location = location;
        this.biopsyLocations = biopsyLocations;
        this.entryStage = entryStage;
    }

    String location() {
        return location;
    }

    String entryStage() {
        return entryStage;
    }

    List<String> biopsyLocations() {
        return biopsyLocations;
    }
}
