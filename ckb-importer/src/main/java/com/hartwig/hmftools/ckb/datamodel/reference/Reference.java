package com.hartwig.hmftools.ckb.datamodel.reference;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.DrugInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.GeneInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Reference {

    public abstract int id();

    @Nullable
    public abstract String pubMedId();

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
    public abstract List<DrugInfo> drug();

    @NotNull
    public abstract List<GeneInfo> gene();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapy();

    @NotNull
    public abstract List<TreatmentApproachInfo> treatmentApproach();

    @NotNull
    public abstract List<VariantInfo> variant();
}
