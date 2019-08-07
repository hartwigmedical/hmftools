package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatch implements KbSpecificObject {

    @NotNull
    public abstract List<MolecularMatchCriteriaUnmet> criteriaUnmet();

    @NotNull
    public abstract List<MolecularMatchPrevalence> prevalence();

    @NotNull
    public abstract String score();

    @NotNull
    public abstract String autoGenerateNarrative();

    @NotNull
    public abstract List<MolecularMatchMutations> mutations();

    @NotNull
    public abstract List<MolecularMatchSource> sources();

    @NotNull
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract List<String> includeCondition0();

    @NotNull
    public abstract List<String> includeCondition1();

    @NotNull
    public abstract String uniqueKey();

    @NotNull
    public abstract String civicValue();

    @NotNull
    public abstract String hashKey();

    @NotNull
    public abstract String regulatoryBodyApproved();

    @NotNull
    public abstract String version();

    @NotNull
    public abstract List<String> includeMutation1();

    @Nullable
    public abstract List<String> includeMutation0();

    @Nullable
    public abstract String guidelineBody();

    @NotNull
    public abstract String regulatoryBody();

    @NotNull
    public abstract String customer();

    @NotNull
    public abstract String direction();

    @NotNull
    public abstract String ampcap();

    @Nullable
    public abstract MolecularMatchAst asts();

    @NotNull
    public abstract List<MolecularMatchVariantInfo> variantInfo();

    @Nullable
    public abstract String guidelineVersion();

    @Nullable
    public abstract String institution();

    @NotNull
    public abstract String tier();

    @NotNull
    public abstract List<MolecularMatchTierExplanation> tierExplanation();

    @NotNull
    public abstract String mvld();

    @NotNull
    public abstract List<MolecularMatchTags> tags();

    @NotNull
    public abstract List<String> criteriaMet();

    @NotNull
    public abstract String biomarkerClass();

    @NotNull
    public abstract List<MolecularMatchClassification> classification();

    @Nullable
    public abstract List<String> includeDrug1();

    @NotNull
    public abstract List<MolecularMatchTherapeuticContext> therapeuticContext();

    @NotNull
    public abstract String sixtier();

    @Nullable
    public abstract String noTherapyAvailable();

    @Nullable
    public abstract String external_id();

    @NotNull
    public abstract String narrative();

    @NotNull
    public abstract String expression();

    @Nullable
    public abstract List<String> includeGene0();

}
