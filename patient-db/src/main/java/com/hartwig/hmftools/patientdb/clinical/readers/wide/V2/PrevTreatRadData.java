package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PrevTreatRadData
{
    @NotNull
    public abstract String subjectKey();
    @NotNull
    public abstract String itemGroupOid();
    @NotNull
    public abstract String itemGroupRepeatKey();
    public abstract String radioSite();
    public abstract String medicalHistoryCategory();
    public abstract int cumulativeDose();
}
