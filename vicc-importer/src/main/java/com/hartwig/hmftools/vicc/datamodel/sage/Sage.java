package com.hartwig.hmftools.vicc.datamodel.sage;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Sage implements KbSpecificObject {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract String clinicalManifestation();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract String evidenceLabel();

    @NotNull
    public abstract String drugLabels();

    @NotNull
    public abstract String germlineOrSomatic();

    @NotNull
    public abstract String publicationUrl();

}
