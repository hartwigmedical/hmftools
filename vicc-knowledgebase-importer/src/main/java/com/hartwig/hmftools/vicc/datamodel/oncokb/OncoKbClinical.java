package com.hartwig.hmftools.vicc.datamodel.oncokb;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncoKbClinical implements KbSpecificObject {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String entrezGeneID();

    @NotNull
    public abstract String isoform();

    @NotNull
    public abstract String refSeq();

    @NotNull
    public abstract OncokbVariant oncokbVariant();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String drugPmids();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String levelLabel();

    @NotNull
    public abstract List<OncoKbDrugAbstracts> oncoKbDrugAbstracts();

}
