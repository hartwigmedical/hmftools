package com.hartwig.hmftools.ckb.reference;

import java.util.List;

import com.hartwig.hmftools.ckb.common.GeneInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.ckb.common.TreatmentApproach;
import com.hartwig.hmftools.ckb.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Reference {

    @NotNull
    public abstract String id();

    @Nullable
    public abstract String pubmedId();

    @Nullable
    public abstract String title();

    @Nullable
    public abstract String url();

    @Nullable
    public abstract String authors();

    @Nullable
    public abstract String journal();

    @Nullable
    public abstract String volume();

    @Nullable
    public abstract String issue();

    @Nullable
    public abstract String date();

    @Nullable
    public abstract String abstractText();

    @Nullable
    public abstract String year();

    @NotNull
    public abstract List<ReferenceDrug> drug();

    @NotNull
    public abstract List<GeneInfo> gene();

    @NotNull
    public abstract List<ReferenceEvidence> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapy();

    @NotNull
    public abstract List<TreatmentApproach> treatmentApproach();

    @NotNull
    public abstract List<VariantInfo> variant();








}
