package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorMeasureData
{
    @NotNull
    public abstract String combinedKey();

    @NotNull
    public abstract String subjectKey();

    @NotNull
    public abstract String formRepeatKey();

    public abstract LocalDate measureDate();

    public abstract String TMDUMTXT(); // TODO find out what this is

    public abstract int recist(); // TODO what is this, very cryptic
    public abstract boolean continueTreatment();

    public abstract EndOfTreatmentReason reasonEndOfTreatment();

    public abstract String reasonEndTreatmentSpecification();


    public enum EndOfTreatmentReason {
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
