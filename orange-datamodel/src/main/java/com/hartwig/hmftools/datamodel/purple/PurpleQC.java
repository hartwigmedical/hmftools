package com.hartwig.hmftools.datamodel.purple;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Set;

@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PurpleQC {
    public abstract Set<PurpleQCStatus> status();
}
