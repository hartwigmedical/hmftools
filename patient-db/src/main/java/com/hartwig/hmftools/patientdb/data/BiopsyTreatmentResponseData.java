package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentResponseData implements Comparable<BiopsyTreatmentResponseData> {

    @Nullable
    public abstract Integer treatmentId();

    @Nullable
    public abstract LocalDate assessmentDate();

    @Nullable
    public abstract LocalDate responseDate();

    @Nullable
    public abstract String response();

    @Nullable
    public abstract String measurementDone();

    @Nullable
    public abstract String boneOnlyDisease();

    @NotNull
    public abstract FormStatus formStatus();

    @NotNull
    public static BiopsyTreatmentResponseData of(@Nullable final LocalDate assessmentDate, @Nullable final LocalDate responseDate,
            @Nullable final String response, @Nullable final String measurementDone, @Nullable final String boneOnlyDisease,
            @NotNull final FormStatus formStatus) {
        return ImmutableBiopsyTreatmentResponseData.of(null,
                assessmentDate,
                responseDate,
                response,
                measurementDone,
                boneOnlyDisease, formStatus);
    }

    @Nullable
    public LocalDate date() {
        if (assessmentDate() == null) {
            return responseDate();
        }
        return assessmentDate();
    }

    public boolean isNotDoneResponse() {
        String response = response();
        return response != null && response.equalsIgnoreCase("nd");
    }

    @Override
    public String toString() {
        return response() + " (" + date() + ")";
    }

    @Override
    public int compareTo(@NotNull final BiopsyTreatmentResponseData other) {
        LocalDate date1 = date();
        LocalDate date2 = other.date();
        if (date1 == null && date2 == null) {
            return 0;
        } else if (date1 == null) {
            return 1;
        } else if (date2 == null) {
            return -1;
        } else {
            return date1.compareTo(date2);
        }
    }
}
