package com.hartwig.hmftools.serve.transvar.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class  TransvarObject {

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String chromosome();

    public abstract long gdnaPosition();

    @NotNull
    public abstract TransvarAnnotation annotation();
}
