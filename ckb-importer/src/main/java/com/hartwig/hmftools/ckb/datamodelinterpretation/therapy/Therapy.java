package com.hartwig.hmftools.ckb.datamodelinterpretation.therapy;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Therapy {

    public abstract int id();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract String synonyms();

    @NotNull
    public abstract List<TherapyDescription> descriptions();

    @NotNull
    public abstract Date createDate();

    @NotNull
    public abstract Date updateDate();

}
