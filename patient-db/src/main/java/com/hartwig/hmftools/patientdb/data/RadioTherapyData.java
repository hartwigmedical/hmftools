package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class RadioTherapyData {
    @Nullable
    private final LocalDate endDate;
    @Nullable
    private final String site;

    public RadioTherapyData(@Nullable final LocalDate endDate, @Nullable final String site) {
        this.endDate = endDate;
        this.site = site;
    }

    public LocalDate endDate() {
        return endDate;
    }

    public String site() {
        return site;
    }
}