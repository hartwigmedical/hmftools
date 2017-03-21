package com.hartwig.hmftools.patientdb;

import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorData {
    private final String location;
    private final List<String> biopsyLocations;
    private final String entryStage;

    TumorData(@Nullable String location, @NotNull List<String> biopsyLocations, @Nullable String entryStage) {
        this.location = location;
        this.biopsyLocations = biopsyLocations;
        this.entryStage = entryStage;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(location).append(" - ").append(entryStage).append(": ").append(biopsyLocations).append("\n");
        return bf.toString();
    }
}
