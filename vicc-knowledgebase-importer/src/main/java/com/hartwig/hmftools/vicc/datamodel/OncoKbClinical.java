package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncoKbClinical implements KbSpecificObject {

    @NotNull
    public abstract String RefSeq();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String Isoform();

    @NotNull
    public abstract OncokbVariant oncokbVariant();

    @NotNull
    public abstract String entrezGeneID();

    @NotNull
    public abstract String drugPmids();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String levelLabel();

    @NotNull
    public abstract List<OncoKbDrugAbstracts> oncoKbDrugAbstracts();


}
