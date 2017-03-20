package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.NotNull;

public class SystemicTherapyData {
    private final Date date;
    private final String type;
    private final String reg;
    private final String bestResponse;

    public SystemicTherapyData(@NotNull Date date, @NotNull String type, @NotNull String reg,
            @NotNull String bestResponse) {
        this.date = date;
        this.type = type;
        this.reg = reg;
        this.bestResponse = bestResponse;
    }
}
