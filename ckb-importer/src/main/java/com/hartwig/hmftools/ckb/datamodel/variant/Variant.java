package com.hartwig.hmftools.ckb.datamodel.variant;

import java.time.LocalDate;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Variant {

    public abstract int id();

    @Nullable
    public abstract Gene gene();

    @NotNull
    public abstract String fullName();

    @Nullable
    public abstract String impact();

    @Nullable
    public abstract String proteinEffect();

    @NotNull
    public abstract List<VariantDescription> variantDescriptions();

    @Nullable
    public abstract String type();

    @NotNull
    public abstract String variant();

    @Nullable
    public abstract LocalDate createDate();

    @Nullable
    public abstract LocalDate updateDate();

    @Nullable
    public abstract ReferenceTranscriptCoordinate referenceTranscriptCoordinate();

    @NotNull
    public abstract List<CategoryVariantPath> categoryVariantPaths();

    @NotNull
    public abstract List<ReferenceTranscriptCoordinate> allTranscriptCoordinates();

    @NotNull
    public abstract List<MemberVariant> memberVariants();
}