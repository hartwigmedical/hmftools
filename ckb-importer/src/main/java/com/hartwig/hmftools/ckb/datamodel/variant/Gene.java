package com.hartwig.hmftools.ckb.datamodel.variant;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Gene {

    public abstract int id();

    @NotNull
    public abstract LocalDate createDate();

    @Nullable
    public abstract LocalDate updateDate();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract String geneRole();

    @Nullable
    public abstract String entrezId();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String mapLocation();

    @Nullable
    public abstract String canonicalTranscript();

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<Reference> references();

}