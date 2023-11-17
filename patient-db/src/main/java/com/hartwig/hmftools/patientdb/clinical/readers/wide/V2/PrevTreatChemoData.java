package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.util.Optional;

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

    @NotNull
    public abstract String itemGroupOid();

    @NotNull
    public abstract Integer itemGroupRepeatedKey();

    public abstract Optional<String> medicalHistoryCategory();

    public abstract Optional<Integer> chemoDrugOne();

    public abstract Optional<Integer> CHTS();

    public abstract Optional<Integer> chemoDrugTwo();

    public abstract Optional<Integer> CHTS1();

    public abstract Optional<Integer> chemoDrugThree();

    public abstract Optional<Integer> CHTS2();

    public abstract Optional<Integer> chemoDrugFour();

    public abstract Optional<Integer> CHTS3();

    public static ImmutablePrevTreatChemoData.Builder builder()
    {
        return ImmutablePrevTreatChemoData.builder();
    }
}
