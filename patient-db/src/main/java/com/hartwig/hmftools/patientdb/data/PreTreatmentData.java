package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class PreTreatmentData implements TreatmentData {

    @Nullable
    public abstract String treatmentGiven();

    @NotNull
    public abstract List<DrugData> drugs();

    @NotNull
    public abstract FormStatusState formStatus();

    public abstract boolean formLocked();

    @Override
    public String toString() {
        return treatmentName() + " (" + startDate() + " - " + endDate() + ")";
    }
}
