package com.hartwig.hmftools.virusinterpreter.algo;

import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxidType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusBlacklistingDb
{

    public abstract int taxid();

    @NotNull
    public abstract TaxidType type();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String reason();
}
