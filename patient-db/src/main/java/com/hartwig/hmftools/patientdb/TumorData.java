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

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(location).append(" - ").append(entryStage).append(": ").append(biopsyLocations).append("\n");
        return bf.toString();
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
