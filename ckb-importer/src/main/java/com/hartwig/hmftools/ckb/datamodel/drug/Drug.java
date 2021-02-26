package com.hartwig.hmftools.ckb.datamodel.drug;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drug {

    public abstract int id();

    @NotNull
    public abstract LocalDate createDate();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<DrugClass> drugClasses();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String tradeName();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String ncitId();

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<Reference> references();
}
