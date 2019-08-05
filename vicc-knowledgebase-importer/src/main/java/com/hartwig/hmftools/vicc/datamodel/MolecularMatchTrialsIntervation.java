package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsIntervation {

    @Nullable
    public abstract String intervention_name();

    @Nullable
    public abstract List<String> other_name();

    @Nullable
    public abstract String description();

    @Nullable
    public abstract List<String> arm_group_label();

    @Nullable
    public abstract String intervention_type();
}
