package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.Nullable;

class RadioTherapyData {
    private final Date endDate;
    private final String site;

    RadioTherapyData(@Nullable Date endDate, @Nullable String site) {
        this.endDate = endDate;
        this.site = site;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(endDate).append(" - ").append(site).append("\n");
        return bf.toString();
    }
}