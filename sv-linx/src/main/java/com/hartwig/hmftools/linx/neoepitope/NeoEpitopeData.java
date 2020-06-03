package com.hartwig.hmftools.linx.neoepitope;

import com.hartwig.hmftools.linx.fusion.GeneFusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class NeoEpitopeData
{
    @NotNull
    public abstract GeneFusion fusion();

    public abstract String upstreamAcids();

    public abstract String downstreamAcids();

    public abstract String novelAcid();

    public abstract int downstreamNmdBases();

}
