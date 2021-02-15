package com.hartwig.hmftools.ckb.datamodel.knowngenomicalteration;

import java.util.List;

import com.hartwig.hmftools.ckb.json.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownGenomicAlteration {

    @NotNull
    public abstract List<Variant> variants();
}