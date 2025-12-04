package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.Driver;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Disruption extends Driver {
    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    String gene();

    @NotNull
    String transcript();

    @NotNull LinxBreakendType type();

    @Nullable
    Double disruptedCopies();

    @Nullable
    Double undisruptedCopies();

    @Nullable
    Integer clusterId();

    @Nullable
    LinxBreakend breakendStart();

    @Nullable
    LinxBreakend breakendEnd();
}
