package com.hartwig.hmftools.ckb.datamodelinterpretation.drug;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drug {

    public abstract int id();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract String tradeName();

    @NotNull
    public abstract List<DrugDescription> drugDescriptions();

    @NotNull
    public abstract String casRegistryNum();

    @NotNull
    public abstract String ncitId();

    @NotNull
    public abstract Date createDate();

}
