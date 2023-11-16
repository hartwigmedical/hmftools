package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TreatChemoAvlData
{
    @NotNull
    abstract String SubjectKey();

    abstract String chemoCode();

    abstract LocalDate TRASDT(); //TODO find better name

    abstract LocalDate TRAEDT(); // TODO find better name
}
