package com.hartwig.hmftools.vicc.datamodel.pmkb;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Pmkb implements KbSpecificObject {

    @NotNull
    public abstract PmkbTumor tumor();

    @NotNull
    public abstract List<PmkbTissue> tissues();

    @NotNull
    public abstract PmkbVariant variant();

}
