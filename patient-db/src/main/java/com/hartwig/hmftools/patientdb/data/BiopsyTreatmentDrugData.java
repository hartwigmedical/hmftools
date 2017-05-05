package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class BiopsyTreatmentDrugData {
    @Nullable
    private final String name;
    @Nullable
    private final String type;
    @Nullable
    private final LocalDate startDate;
    @Nullable
    private final LocalDate endDate;

    public BiopsyTreatmentDrugData(@Nullable final String name, @Nullable final String type,
            @Nullable final LocalDate startDate, @Nullable final LocalDate endDate) {
        this.name = name;
        this.type = type;
        this.startDate = startDate;
        this.endDate = endDate;
    }

    @Nullable
    public String name() {
        return name;
    }

    @Nullable
    public String type() {
        return type;
    }

    @Nullable
    public LocalDate startDate() {
        return startDate;
    }

    @Nullable
    public LocalDate endDate() {
        return endDate;
    }
}
