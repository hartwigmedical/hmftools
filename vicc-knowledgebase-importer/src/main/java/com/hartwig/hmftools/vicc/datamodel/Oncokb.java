package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Oncokb implements KbSpecificObject {

    @Nullable
    public abstract BiologicalOncoKb biologicalOncoKb();

    @Nullable
    public abstract ClinicalOncoKb clinicalOncoKb();

}
