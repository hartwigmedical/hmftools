package com.hartwig.hmftools.ckb.datamodel.drug;

import java.time.LocalDate;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drug {

    public abstract int id();

    @Nullable
    public abstract LocalDate createDate();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<DrugClass> drugClasses(); //some drugs has none drug class eg. drugs 9758

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DrugDescription> descriptions();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String ncitId();
}
