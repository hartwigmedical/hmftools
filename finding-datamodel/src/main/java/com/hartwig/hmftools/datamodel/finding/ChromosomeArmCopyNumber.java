package com.hartwig.hmftools.datamodel.finding;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChromosomeArmCopyNumber extends Finding
{
    enum ChromosomeArm
    {
        P_ARM,
        Q_ARM,
    }

    @NotNull
    String chromosome();

    @NotNull
    ChromosomeArm arm();

    double meanCopyNumber();
}
