package com.hartwig.hmftools.ckb.datamodelinterpretation.variant;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Variant {

    public abstract int id();

    @NotNull
    public abstract String fullName();

    @Nullable
    public abstract String impact();

    @Nullable
    public abstract String proteinEffect();

//    @NotNull
//    public abstract List<VariantDescription> variantDescriptions();
//
//    @NotNull
//    public abstract String type();
//
//    @NotNull
//    public abstract String variant();
//
//    @NotNull
//    public abstract Date createDate();
//
//    @NotNull
//    public abstract Date updateDate();
//
//    @NotNull
//    public abstract ReferenceTranscriptCoordinate referenceTranscriptCoordinate();
//
//    @NotNull
//    public abstract List<CategoryVariantPath> categoryVariantPaths();
//
//    @NotNull
//    public abstract List<ReferenceTranscriptCoordinate> allTranscriptCoordinated();
//
//    @NotNull
//    public abstract List<MemberVariant> memberVariants();

}
