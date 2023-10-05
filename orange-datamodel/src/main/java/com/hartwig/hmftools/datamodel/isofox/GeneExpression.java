package com.hartwig.hmftools.datamodel.isofox;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneExpression
{
    @NotNull
    String gene();

    double tpm();

    double medianTpmCancer();

    double percentileCancer();

    @Nullable
    Double medianTpmCohort();

    @Nullable
    Double percentileCohort();
}
