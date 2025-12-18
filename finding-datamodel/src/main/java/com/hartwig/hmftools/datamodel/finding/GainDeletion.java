package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GainDeletion extends Driver
{
    enum Type
    {
        GERMLINE_DEL_HOM_IN_TUMOR,
        GERMLINE_DEL_HET_IN_TUMOR,
        SOMATIC_GAIN,
        SOMATIC_DEL
    }

    @NotNull
    Type type();

    @NotNull
    CopyNumberInterpretation interpretation();

    @NotNull
    String gene();

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    String transcript();

    boolean isCanonical();

    double minCopies();

    double maxCopies();
}
