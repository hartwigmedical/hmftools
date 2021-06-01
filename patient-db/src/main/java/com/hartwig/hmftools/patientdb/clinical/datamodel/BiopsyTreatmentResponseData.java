package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

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
    public static BiopsyTreatmentResponseData of(@Nullable Integer treatmentId, @Nullable LocalDate assessmentDate,
            @Nullable LocalDate responseDate, @Nullable String response, @Nullable String measurementDone, @Nullable String boneOnlyDisease,
            @NotNull FormStatus formStatus) {
        return ImmutableBiopsyTreatmentResponseData.builder()
                .treatmentId(treatmentId)
                .assessmentDate(assessmentDate)
                .responseDate(responseDate)
                .response(response)
                .measurementDone(measurementDone)
                .boneOnlyDisease(boneOnlyDisease)
                .formStatus(formStatus)
                .build();
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
    public int compareTo(@NotNull BiopsyTreatmentResponseData other) {
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
