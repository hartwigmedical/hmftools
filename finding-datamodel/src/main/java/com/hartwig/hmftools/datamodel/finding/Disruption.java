package com.hartwig.hmftools.datamodel.finding;

import java.util.Set;

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

        private static final Set<Type> SOMATIC_TYPES = Set.of(SOMATIC_DISRUPTION, SOMATIC_HOM_DUP_DISRUPTION, SOMATIC_HOM_DEL_DISRUPTION);
        private static final Set<Type> HOMOZYGOUS_TYPES = Set.of(SOMATIC_HOM_DUP_DISRUPTION, SOMATIC_HOM_DEL_DISRUPTION, GERMLINE_HOM_DUP_DISRUPTION);

        public boolean isSomatic() { return SOMATIC_TYPES.contains(this); }
        public boolean isGermline() { return !SOMATIC_TYPES.contains(this); }
        public boolean isHomozygous() { return HOMOZYGOUS_TYPES.contains(this); }
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
