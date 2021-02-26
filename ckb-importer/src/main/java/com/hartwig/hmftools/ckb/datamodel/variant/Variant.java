package com.hartwig.hmftools.ckb.datamodel.variant;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Variant {

    public abstract int id();

    @NotNull
    public abstract LocalDate createDate();

    @NotNull
    public abstract LocalDate updateDate();

    @NotNull
    public abstract String fullName();

    @NotNull
    public abstract String variant();

    @Nullable
    public abstract String impact();

    @Nullable
    public abstract String proteinEffect();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract Gene gene();

    @Nullable
    public abstract TranscriptCoordinate referenceTranscriptCoordinate();

    @NotNull
    public abstract List<TranscriptCoordinate> allTranscriptCoordinates();

    @NotNull
    public abstract List<String> categoryVariantPaths();

    @NotNull
    public abstract List<MemberVariant> memberVariants();

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<Reference> references();

}