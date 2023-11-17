package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;
import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class })
public abstract class TreatChemoAvlData
{
    @NotNull
    abstract String subjectKey();

    abstract Optional<String> chemoCode();

    abstract Optional<LocalDate> TRASDT(); //TODO find better name

    abstract Optional<LocalDate> TRAEDT(); // TODO find better name

    public static ImmutableTreatChemoAvlData.Builder builder() {
        return ImmutableTreatChemoAvlData.builder();
    }
}
