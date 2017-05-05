package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.base.Strings;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyTreatmentData {

    private final int id;
    @Nullable
    private final String treatmentGiven;
    @Nullable
    private final LocalDate startDate;
    @Nullable
    private final LocalDate endDate;
    @NotNull
    private final List<BiopsyTreatmentDrugData> drugs;
    @Nullable
    private final Integer biopsyId;

    private static final AtomicInteger idCounter = new AtomicInteger();

    private static int createId() {
        return idCounter.getAndIncrement();
    }

    public BiopsyTreatmentData(final int id, @Nullable final String treatmentGiven,
            @Nullable final LocalDate startDate, @Nullable final LocalDate endDate,
            @NotNull final List<BiopsyTreatmentDrugData> drugs, @Nullable final Integer biopsyId) {
        this.id = id;
        this.treatmentGiven = treatmentGiven;
        this.startDate = startDate;
        this.endDate = endDate;
        this.drugs = drugs;
        this.biopsyId = biopsyId;
    }

    public BiopsyTreatmentData(@Nullable final String treatmentGiven, @Nullable final LocalDate startDate,
            @Nullable final LocalDate endDate, @NotNull final List<BiopsyTreatmentDrugData> drugs) {
        this(createId(), treatmentGiven, startDate, endDate, drugs, null);
    }

    @Nullable
    public LocalDate startDate() {
        return startDate;
    }

    @Nullable
    public LocalDate endDate() {
        return endDate;
    }

    @NotNull
    public List<BiopsyTreatmentDrugData> drugs() {
        return drugs;
    }

    @Nullable
    public String treatmentName() {
        final StringJoiner joiner = new StringJoiner("/");
        drugs.forEach(drug -> {
            final String drugName = drug.name();
            if (drugName == null) {
                joiner.add("NULL");
            } else {
                joiner.add(drugName);
            }
        });
        return Strings.emptyToNull(joiner.toString());
    }

    @Nullable
    public String type() {
        final StringJoiner joiner = new StringJoiner("/");
        drugs.forEach(drug -> {
            final String drugType = drug.type();
            if (drugType == null) {
                joiner.add("NULL");
            } else {
                joiner.add(drugType);
            }
        });
        return Strings.emptyToNull(joiner.toString());
    }

    public int id() {
        return id;
    }

    @Nullable
    public Integer biopsyId() {
        return biopsyId;
    }

    @Nullable
    public String treatmentGiven() {
        return treatmentGiven;
    }
}
