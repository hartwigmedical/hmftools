package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Pmkb implements KbSpecificObject {

    @NotNull
    public abstract List<TumorPmkb> tumor();

    @NotNull
    public abstract List<TissuePmkb> tissue();

    @NotNull
    public abstract List<VariantPmkb> variant();


}
