package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

class RadioTherapyData {
    private final LocalDate endDate;
    private final String site;

    RadioTherapyData(@Nullable LocalDate endDate, @Nullable String site) {
        this.endDate = endDate;
        this.site = site;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(endDate).append(" - ").append(site).append("\n");
        return bf.toString();
    }

    LocalDate endDate() {
        return endDate;
    }

    String site() {
        return site;
    }
}