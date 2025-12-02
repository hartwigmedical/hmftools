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
public interface PurpleDriver
{
    @NotNull
    String gene();

    @NotNull
    String transcript();

    @NotNull
    PurpleDriverType type();

    double driverLikelihood();

    @NotNull
    PurpleLikelihoodMethod likelihoodMethod();

    @NotNull
    ReportedStatus reportedStatus();

    @NotNull
    DriverInterpretation driverInterpretation();

    boolean isCanonical();
}
