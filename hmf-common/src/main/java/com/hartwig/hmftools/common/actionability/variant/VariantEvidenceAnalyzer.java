package com.hartwig.hmftools.common.actionability.variant;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VariantEvidenceAnalyzer {

    private static final Set<CodingEffect> CODING_EFFECTS =
            Sets.newHashSet(CodingEffect.SPLICE, CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE);

    @NotNull
    private final List<ActionableVariant> actionableVariants;
    @NotNull
    private final List<ActionableRange> actionableRanges;

    VariantEvidenceAnalyzer(@NotNull final List<ActionableVariant> actionableVariants, @NotNull List<ActionableRange> actionableRanges) {
        this.actionableVariants = actionableVariants;
        this.actionableRanges = actionableRanges;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        actionableRanges.forEach(x -> genes.add(x.gene()));
        actionableVariants.forEach(x -> genes.add(x.gene()));
        return genes;
    }

    @NotNull
    public List<ActionableVariant> actionableVariants() {
        return actionableVariants;
    }

    @NotNull
    public List<ActionableRange> actionableRanges() {
        return actionableRanges;
    }

    @NotNull
    public List<EvidenceItem> evidenceForVariant(@NotNull Variant variant, @Nullable String primaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (ActionableVariant actionableVariant : actionableVariants) {
            if (variant.gene().equals(actionableVariant.gene()) && variant.chromosome().equals(actionableVariant.chromosome())
                    && variant.position() == actionableVariant.position() && variant.ref().equals(actionableVariant.ref()) && variant.alt()
                    .equals(actionableVariant.alt())) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableVariant(actionableVariant);
                evidenceBuilder.event(eventString(variant));
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.isCancerTypeMatch(actionableVariant.cancerType(), primaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }

        for (ActionableRange actionableRange : actionableRanges) {
            if (CODING_EFFECTS.contains(variant.canonicalCodingEffect()) && variant.gene().equals(actionableRange.gene())
                    && variant.chromosome().equals(actionableRange.chromosome()) && variant.position() >= actionableRange.start()
                    && variant.position() <= actionableRange.end()) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableRange(actionableRange);
                evidenceBuilder.event(eventString(variant));
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.isCancerTypeMatch(actionableRange.cancerType(), primaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }
        return evidenceItems;
    }

    @NotNull
    public static String eventString(@NotNull Variant variant) {
        String description = variant.canonicalCodingEffect() == CodingEffect.SPLICE
                ? variant.canonicalHgvsCodingImpact()
                : variant.canonicalHgvsProteinImpact();
        return variant.gene() + " " + description;
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableVariant(@NotNull ActionableVariant actionableVariant) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableVariant.reference())
                .source(ActionabilitySource.fromString(actionableVariant.source()))
                .drug(actionableVariant.drug())
                .drugsType(actionableVariant.drugsType())
                .level(EvidenceLevel.fromString(actionableVariant.level()))
                .response(actionableVariant.response())
                .cancerType(actionableVariant.cancerType())
                .scope(EvidenceScope.SPECIFIC);
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableRange(@NotNull ActionableRange actionableRange) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableRange.reference())
                .source(ActionabilitySource.fromString(actionableRange.source()))
                .drug(actionableRange.drug())
                .drugsType(actionableRange.drugsType())
                .level(EvidenceLevel.fromString(actionableRange.level()))
                .response(actionableRange.response())
                .cancerType(actionableRange.cancerType())
                .scope(EvidenceScope.GENE_LEVEL);
    }
}