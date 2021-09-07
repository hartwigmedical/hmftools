package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.util.List;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PreTreatmentData implements TreatmentData {

    @Nullable
    @Override
    public abstract String treatmentGiven();

    @Nullable
    @Override
    public abstract String radiotherapyGiven();

    @NotNull
    @Override
    public abstract List<DrugData> drugs();

    @NotNull
    @Override
    public abstract FormStatus formStatus();

    @Override
    public String toString() {
        return treatmentName() + " (" + startDate() + " - " + endDate() + ")";
    }
}
