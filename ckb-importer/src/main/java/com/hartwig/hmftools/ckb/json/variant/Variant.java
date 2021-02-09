package com.hartwig.hmftools.ckb.json.variant;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GeneInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Variant implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String fullName();

    @Nullable
    public abstract String impact();

    @Nullable
    public abstract String proteinEffect();

    @NotNull
    public abstract List<DescriptionInfo> descriptions();

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
    public abstract List<VariantPartnerGene> partnerGenes();

    @NotNull
    public abstract List<VariantCategoryVariantPath> categoryVariantPaths();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<EvidenceInfo> extendedEvidence();

    @NotNull
    public abstract List<MolecularProfileInfo> molecularProfiles();

    @NotNull
    public abstract List<VariantTranscriptCoordinate> allTranscriptCoordinates();

    @NotNull
    public abstract List<VariantInfo> memberVariants();
}
