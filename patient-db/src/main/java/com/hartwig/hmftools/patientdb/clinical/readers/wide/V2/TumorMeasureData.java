package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;
import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class })
public abstract class TumorMeasureData
{
    @NotNull
    public abstract String combinedKey();

    @NotNull
    public abstract String subjectKey();

    @NotNull
    public abstract String formRepeatKey();

    public abstract Optional<LocalDate> measureDate();

    public abstract Optional<String> TMDUMTXT(); // TODO find out what this is

    public abstract Optional<Integer> recist(); // TODO what is this, very cryptic

    public abstract Optional<ResponseEvaluation> responseEvaluation();

    public abstract Optional<EndOfTreatmentReason> reasonEndOfTreatment();

    public abstract Optional<String> reasonEndOfTreatmentSpecification();

    public static ImmutableTumorMeasureData.Builder builder() {
        return ImmutableTumorMeasureData.builder();
    }

    public enum ResponseEvaluation {
        CLINICAL_BENEFIT,
        STOP_TREATMENT
    }

    public enum EndOfTreatmentReason
    {
        CLINICAL_DETERIORATION,
        DEATH,
        RADIOLOGIC_PROGRESSION_NON_RECIST_LESIONS,
        VISIBLE_EVIDENT_PROGRESSION,
        TUMOR_MARKERS,
        CHEMO_PAUSE,
        CONFORM_PROTOCOL_CURATIVE_SETTING,
        OTHER
    }

}
