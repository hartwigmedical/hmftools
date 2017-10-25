package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData {

    public abstract int id();

    @Nullable
    public abstract String treatmentGiven();

    @Nullable
    public abstract LocalDate startDate();

    @Nullable
    public abstract LocalDate endDate();

    @NotNull
    public abstract List<BiopsyTreatmentDrugData> drugs();

    @Nullable
    public abstract Integer biopsyId();

    @NotNull
    public abstract String formStatus();

    @NotNull
    public abstract String formLocked();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable final String treatmentGiven, @Nullable final LocalDate startDate,
            @Nullable final LocalDate endDate, @NotNull final List<BiopsyTreatmentDrugData> drugs, @NotNull final String formStatus,
            @NotNull final String formLocked) {

        return ImmutableBiopsyTreatmentData.of(createId(), treatmentGiven, startDate, endDate, drugs, null, formStatus, formLocked);
    }

    @Nullable
    public String treatmentName() {
        final List<String> drugNames = Lists.newArrayList();
        for (BiopsyTreatmentDrugData drug : drugs()) {
            final String drugName = drug.name();
            if (drugName == null) {
                drugNames.add("NULL");
            } else {
                drugNames.add(toFirstLetterUpperCase(drugName));
            }
        }

        Collections.sort(drugNames);

        final StringJoiner joiner = new StringJoiner("/");
        for (String drugName : drugNames) {
            joiner.add(drugName);
        }
        return Strings.emptyToNull(joiner.toString());
    }

    @Nullable
    public String type() {
        final StringJoiner joiner = new StringJoiner("/");
        drugs().forEach(drug -> {
            final String drugType = drug.type();
            if (drugType == null) {
                joiner.add("NULL");
            } else {
                joiner.add(drugType);
            }
        });
        return Strings.emptyToNull(joiner.toString());
    }

    @Override
    public String toString() {
        return treatmentName() + "(" + startDate() + " - " + endDate() + ")";
    }

    @NotNull
    private static String toFirstLetterUpperCase(@NotNull final String string) {
        if (string.length() == 0) {
            return string;
        } else if (string.length() == 1) {
            return string.toUpperCase();
        } else {
            return string.toUpperCase().substring(0, 1) + string.substring(1);
        }
    }
}
