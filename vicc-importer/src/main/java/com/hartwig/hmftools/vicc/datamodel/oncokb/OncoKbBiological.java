package com.hartwig.hmftools.vicc.datamodel.oncokb;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncoKbBiological implements KbSpecificObject {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String entrezGeneId();

    @NotNull
    public abstract String isoform();

    @NotNull
    public abstract String refSeq();

    @NotNull
    public abstract OncoKbVariant oncokbVariant();

    @NotNull
    public abstract String oncogenic();

    @NotNull
    public abstract String mutationEffect();

    @NotNull
    public abstract String mutationEffectPmids();

    @NotNull
    public abstract String mutationEffectAbstracts();

}
