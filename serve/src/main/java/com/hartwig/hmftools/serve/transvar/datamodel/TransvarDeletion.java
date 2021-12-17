package com.hartwig.hmftools.serve.transvar.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TransvarDeletion implements TransvarAnnotation {

    public abstract int deletedBaseCount();

    public abstract int leftAlignedGDNAPosition();
}
