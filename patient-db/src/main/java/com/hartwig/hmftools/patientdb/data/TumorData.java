package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import org.jetbrains.annotations.Nullable;

public class TumorData {
    @Nullable
    private final String location;
    @Nullable
    private final List<String> biopsyLocations;
    @Nullable
    private final String entryStage;

    public TumorData(@Nullable final String location, @Nullable final List<String> biopsyLocations,
            @Nullable final String entryStage) {
        this.location = location;
        this.biopsyLocations = biopsyLocations;
        this.entryStage = entryStage;
    }

    public String location() {
        return location;
    }

    public String entryStage() {
        return entryStage;
    }

    public List<String> biopsyLocations() {
        return biopsyLocations;
    }
}
