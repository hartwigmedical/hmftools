package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatch implements KbSpecificObject {

    @NotNull
    public abstract String direction();

    @NotNull
    public abstract String biomarkerClass();

    @NotNull
    public abstract List<MolecularMatchMutation> mutations();

    @NotNull
    public abstract List<MolecularMatchVariantInfo> variantInfos();

    @NotNull
    public abstract List<MolecularMatchPrevalence> prevalences();

    @NotNull
    public abstract String score();

    @NotNull
    public abstract List<MolecularMatchSource> sources();

    @NotNull
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String tier();

    @NotNull
    public abstract List<MolecularMatchTierExplanation> tierExplanations();

    @NotNull
    public abstract String ampcap();

    @NotNull
    public abstract String civicValue();

    @NotNull
    public abstract String regulatoryBody();

    @NotNull
    public abstract String regulatoryBodyApproved();

    @Nullable
    public abstract String guidelineBody();

    @Nullable
    public abstract String guidelineVersion();

    @NotNull
    public abstract List<String> includeGene1();

    @NotNull
    public abstract List<String> includeFinding1();

    @NotNull
    public abstract List<String> includeCondition1();

    @NotNull
    public abstract List<String> includeMutation1();

    @NotNull
    public abstract List<String> includeDrug1();

    @NotNull
    public abstract List<String> includeDrugClass1();

    @NotNull
    public abstract List<String> includeResistance1();

    @NotNull
    public abstract List<String> includeStage0();

    @NotNull
    public abstract List<String> includeGene0();

    @NotNull
    public abstract List<String> includeCondition0();

    @NotNull
    public abstract List<String> includeMutation0();

    @NotNull
    public abstract List<String> criteriaMets();

    @NotNull
    public abstract List<MolecularMatchCriteriaUnmet> criteriaUnmets();

    @NotNull
    public abstract MolecularMatchAst ast();

    @NotNull
    public abstract List<String> institutions();

    @NotNull
    public abstract List<MolecularMatchTag> tags();

    @NotNull
    public abstract List<MolecularMatchClassification> classifications();

    @Nullable
    public abstract String noTherapyAvailable();

    @NotNull
    public abstract List<MolecularMatchTherapeuticContext> therapeuticContexts();

    @NotNull
    public abstract String sixtier();

    @NotNull
    public abstract String mvld();

    @NotNull
    public abstract String autoGenerateNarrative();

    @NotNull
    public abstract String narrative();

    @NotNull
    public abstract String expression();

    @NotNull
    public abstract String customer();

    @NotNull
    public abstract String version();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract List<String> externalIds();

    @NotNull
    public abstract String uniqueKey();

    @NotNull
    public abstract String hashKey();

}
