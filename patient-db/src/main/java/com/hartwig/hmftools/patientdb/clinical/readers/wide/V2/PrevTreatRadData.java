package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class })
public abstract class PrevTreatRadData
{
    @NotNull
    public abstract String subjectKey();

    @NotNull
    public abstract String itemGroupOid();

    @NotNull
    public abstract Integer itemGroupRepeatKey();

    public abstract Optional<String> radioSite();

    public abstract Optional<String> medicalHistoryCategory();

    public abstract Optional<Integer> cumulativeDose();

    public static ImmutablePrevTreatRadData.Builder builder()
    {
        return ImmutablePrevTreatRadData.builder();
    }
}
