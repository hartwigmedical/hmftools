package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleGermlineDeletion
{
    @NotNull
    String gene();

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    int regionStart();

    int regionEnd();

    int depthWindowCount();

    int exonStart();

    int exonEnd();

    @NotNull
    PurpleGermlineDetectionMethod detectionMethod();

    @NotNull
    PurpleGermlineStatus normalStatus();

    @NotNull
    PurpleGermlineStatus tumorStatus();

    double germlineCopyNumber();

    double tumorCopyNumber();

    @NotNull
    String filter();

    int cohortFrequency();

    @NotNull
    ReportedStatus reportedStatus();

    @NotNull
    DriverInterpretation driverInterpretation();
}
