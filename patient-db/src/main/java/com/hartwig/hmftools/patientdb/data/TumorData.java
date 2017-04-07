package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorData {
    @Nullable
    private final String location;
    @NotNull
    private final List<String> biopsyLocations;
    @Nullable
    private final String entryStage;

    public TumorData(@Nullable final String location, @NotNull final List<String> biopsyLocations,
            @Nullable final String entryStage) {
        this.location = location;
        this.biopsyLocations = biopsyLocations;
        this.entryStage = entryStage;
    }

    @Nullable
    public String location() {
        return location;
    }

    @Nullable
    public String entryStage() {
        return entryStage;
    }

    @NotNull
    public List<String> biopsyLocations() {
        return biopsyLocations;
    }
}
