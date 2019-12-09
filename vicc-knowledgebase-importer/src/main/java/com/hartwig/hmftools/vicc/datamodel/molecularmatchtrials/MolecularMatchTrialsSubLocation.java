package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsSubLocation {

    @NotNull
    public abstract String type();

    @NotNull
    public abstract List<String> coordinates();

}
