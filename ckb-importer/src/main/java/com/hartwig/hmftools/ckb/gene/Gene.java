package com.hartwig.hmftools.ckb.gene;

import java.util.List;

import com.hartwig.hmftools.ckb.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.common.EffectInfo;
import com.hartwig.hmftools.ckb.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Gene {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> term();

    @Nullable
    public abstract String entrezId();

    @NotNull
    public abstract List<String> synonym();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String mapLocation();

    @NotNull
    public abstract List<DescriptionInfo> description();

    @Nullable
    public abstract String canonicalTranscript();

    @NotNull
    public abstract String geneRole();

    @NotNull
    public abstract String createDate();

    @Nullable
    public abstract String updateDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrial();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<VariantInfo> variant();

    @NotNull
    public abstract List<MolecularProfileInfo> molecularProfile();

    @NotNull
    public abstract List<EffectInfo> categoryVariant();


}
