package com.hartwig.hmftools.common.drivercatalog;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverCatalog
{
    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    String gene();

    @NotNull
    DriverCategory category();

    @NotNull
    DriverType driver();

    @NotNull
    LikelihoodMethod likelihoodMethod();

    double driverLikelihood();

    int missense();

    int nonsense();

    int splice();

    int inframe();

    int frameshift();

    boolean biallelic();

    double minCopyNumber();

    double maxCopyNumber();
}
