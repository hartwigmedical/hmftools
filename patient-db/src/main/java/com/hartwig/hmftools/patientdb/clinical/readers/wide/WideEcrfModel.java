package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideEcrfModel {

    @NotNull
    public abstract List<WidePreAvlTreatmentData> preAvlTreatments();

    @NotNull
    public abstract List<WideBiopsyData> biopsies();

    @NotNull
    public abstract List<WideAvlTreatmentData> avlTreatments();

    @NotNull
    public abstract List<WideResponseData> responses();

    @NotNull
    public abstract List<WideFiveDays> fiveDays();
}
