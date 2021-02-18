package com.hartwig.hmftools.ckb.json.gene;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonGene implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> terms();

    @Nullable
    public abstract String entrezId();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String mapLocation();

    @NotNull
    public abstract List<DescriptionInfo> descriptions();

    @Nullable
    public abstract String canonicalTranscript();

    @NotNull
    public abstract String geneRole();

    @NotNull
    public abstract LocalDate createDate();

    @NotNull
    public abstract LocalDate updateDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrials();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<VariantInfo> variants();

    @NotNull
    public abstract List<MolecularProfileInfo> molecularProfiles();

    @NotNull
    public abstract List<VariantInfo> categoryVariants();
}
