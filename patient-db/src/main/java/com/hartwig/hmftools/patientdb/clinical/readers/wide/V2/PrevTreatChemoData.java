package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PrevTreatChemoData // TODO CHTS meaning
{
    @NotNull
    public abstract String subjectKey();
    public abstract String ItemGroupOid();
    public abstract int ItemGroupRepeatedKey();
    public abstract String medicalHistoryCategory();

    public abstract int chemoDrugOne();

    public abstract int CHTS();

    public abstract int chemoDrugTwo();

    public abstract int CHTS1();

    public abstract int chemoDrugThree();

    public abstract int CHTS2();

    public abstract int chemoDrugFour();

    public abstract int CHTS3();
}
