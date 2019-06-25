package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiologicalOncoKb implements KbSpecificObject {

    @NotNull
    public abstract String mutationEffectPmids();

    @NotNull
    public abstract String Isoform();

    @NotNull
    public abstract VariantOncokb variantOncokb();

    @NotNull
    public abstract String entrezGeneID();

    @NotNull
    public abstract String oncogenic();

    @NotNull
    public abstract String mutationEffect();

    @NotNull
    public abstract String RefSeq();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String mutationEffectAbstracts();



}
