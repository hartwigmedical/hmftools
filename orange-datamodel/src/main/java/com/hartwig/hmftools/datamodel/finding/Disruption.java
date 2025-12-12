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

    enum Type
    {
        SOMATIC_DISRUPTION,
        SOMATIC_HOM_DUP_DISRUPTION,
        SOMATIC_HOM_DEL_DISRUPTION,
        GERMLINE_DISRUPTION,
        GERMLINE_HOM_DUP_DISRUPTION;

        boolean isSomatic() { return this == SOMATIC_DISRUPTION || this == SOMATIC_HOM_DUP_DISRUPTION || this == SOMATIC_HOM_DEL_DISRUPTION; }
        boolean isGermline() { return this == GERMLINE_DISRUPTION || this == GERMLINE_HOM_DUP_DISRUPTION; }
        boolean isHomozygous() { return this == SOMATIC_HOM_DEL_DISRUPTION; }
    }

    @NotNull
    Type type();

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    String gene();

    @NotNull
    String transcript();

    boolean isCanonical();

    @NotNull LinxBreakendType breakendType();

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
