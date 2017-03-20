package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.NotNull;

public class RadioTherapyData {
    private final Date date;
    private final String site;

    public RadioTherapyData(@NotNull Date date, @NotNull String site) {
        this.date = date;
        this.site = site;
    }
}