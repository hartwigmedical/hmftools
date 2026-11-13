package com.hartwig.hmftools.datamodel.viridian;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VirusReport
{
    @NotNull
    String name();

    @NotNull
    String virusGroup();

    @NotNull
    ReportQcStatus qcStatus();

    int integrations();

    @Nullable
    VirusInterpretation interpretation();

    double percentageCovered();

    double meanCoverage();

    boolean reported();

    @NotNull
    ClonalityStatus clonalityStatus();

    @Nullable
    Double copiesPerTumorCell();
}
