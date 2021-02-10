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

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DrugDescription> drugDescriptions();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String ncitId();

    @Nullable
    public abstract Date createDate();

}
