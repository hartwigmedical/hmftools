package com.hartwig.hmftools.ckb.datamodel.variant;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EffectInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.GeneInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;

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

    @NotNull
    public abstract List<DescriptionInfo> description();

    @Nullable
    public abstract String type();

    @NotNull
    public abstract GeneInfo gene();

    @NotNull
    public abstract String variant();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();

    @Nullable
    public abstract VariantTranscriptCoordinate referenceTranscriptCoordinate();

    @NotNull
    public abstract List<VariantPartnerGene> partnerGene();

    @NotNull
    public abstract List<VariantCategoryVariantPath> categoryVariantPath();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<EvidenceInfo> extendedEvidence();

    @NotNull
    public abstract List<MolecularProfileInfo> molecularProfile();

    @NotNull
    public abstract List<VariantTranscriptCoordinate> allTranscriptCoordinate();

    @NotNull
    public abstract List<EffectInfo> memberVariant();


}
