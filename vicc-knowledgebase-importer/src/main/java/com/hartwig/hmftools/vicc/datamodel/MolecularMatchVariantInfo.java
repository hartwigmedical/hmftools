package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchVariantInfo {

    @NotNull
    public abstract String classification();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract List<String> consequences();

    @NotNull
    public abstract List<String> fusions();

    @NotNull
    public abstract List<MolecularMatchLocations> locations();

    @NotNull
    public abstract String geneFusionPartner();

    @Nullable
    public abstract String COSMIC_ID();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String popFreqMax();
}
