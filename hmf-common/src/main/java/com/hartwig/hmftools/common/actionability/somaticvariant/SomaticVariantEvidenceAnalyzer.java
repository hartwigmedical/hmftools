package com.hartwig.hmftools.common.actionability.somaticvariant;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticVariantEvidenceAnalyzer {

    private static final Set<CodingEffect> CODING_EFFECTS =
            Sets.newHashSet(CodingEffect.SPLICE, CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE);

    @NotNull
    private final List<ActionableSomaticVariant> actionableVariants;
    @NotNull
    private final List<ActionableRange> actionableRanges;

    SomaticVariantEvidenceAnalyzer(@NotNull final List<ActionableSomaticVariant> actionableVariants,
            @NotNull List<ActionableRange> actionableRanges) {
        this.actionableVariants = actionableVariants;
        this.actionableRanges = actionableRanges;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (ActionableSomaticVariant variant : actionableVariants) {
            genes.add(variant.gene());
        }
        for (ActionableRange range : actionableRanges) {
            genes.add(range.gene());
        }
        return genes;
    }

    @NotNull
    public List<EvidenceItem> evidenceForSomaticVariant(@NotNull SomaticVariant variant, @Nullable String doidsPrimaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (ActionableSomaticVariant actionableVariant : actionableVariants) {
            if (variant.gene().equals(actionableVariant.gene()) && variant.chromosome().equals(actionableVariant.chromosome())
                    && variant.position() == actionableVariant.position() && variant.ref().equals(actionableVariant.ref()) && variant.alt()
                    .equals(actionableVariant.alt())) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableVariant(actionableVariant);

                evidenceBuilder.event(eventString(variant));
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.foundTumorLocation(actionableVariant.cancerType(), doidsPrimaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }

        for (ActionableRange actionableRange : actionableRanges) {
            if (CODING_EFFECTS.contains(variant.canonicalCodingEffect()) && variant.gene().equals(actionableRange.gene())
                    && variant.chromosome().equals(actionableRange.chromosome()) && variant.position() >= actionableRange.start()
                    && variant.position() <= actionableRange.end()) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableRange(actionableRange);

                evidenceBuilder.event(eventString(variant));
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.foundTumorLocation(actionableRange.cancerType(), doidsPrimaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }
        return evidenceItems;
    }

    @NotNull
    private static String eventString(@NotNull SomaticVariant variant) {
        return variant.gene() + " " + variant.canonicalHgvsProteinImpact();
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableVariant(@NotNull ActionableSomaticVariant actionableVariant) {
        return ImmutableEvidenceItem.builder()
                .reference(actionableVariant.reference())
                .source(ActionabilitySource.fromString(actionableVariant.source()))
                .drug(actionableVariant.drug())
                .drugsType(actionableVariant.drugsType())
                .level(EvidenceLevel.fromString(actionableVariant.level()))
                .response(actionableVariant.response())
                .cancerType(actionableVariant.cancerType());
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
                .cancerType(actionableRange.cancerType());
    }
}