package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.confidenceregions.ConfidenceRegionsRule;

import org.jetbrains.annotations.NotNull;

public class VariantAnalyzer {

    @NotNull
    private final ConfidenceRegionsRule confidenceRegionsRule;
    @NotNull
    private final ConsequenceDeterminer determiner;

    public static VariantAnalyzer fromSlicingRegions(@NotNull final GeneModel geneModel, @NotNull final Slicer giabHighConfidenceRegion,
            @NotNull final Slicer cpctSlicingRegion) {
        final ConfidenceRegionsRule confidenceRegionsRule = ConfidenceRegionsRule.fromSlicers(giabHighConfidenceRegion, cpctSlicingRegion);
        final ConsequenceDeterminer determiner = fromHmfSlicingRegion(geneModel);
        return new VariantAnalyzer(confidenceRegionsRule, determiner);
    }

    @NotNull
    private static ConsequenceDeterminer fromHmfSlicingRegion(@NotNull final GeneModel geneModel) {
        return new ConsequenceDeterminer(geneModel);
    }

    private VariantAnalyzer(@NotNull final ConfidenceRegionsRule confidenceRegionsRule, @NotNull final ConsequenceDeterminer determiner) {
        this.confidenceRegionsRule = confidenceRegionsRule;
        this.determiner = determiner;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> passedVariants = passOnly(variants);
        final List<SomaticVariant> consensusPassedVariants = confidenceRegionsRule.removeUnreliableVariants(passedVariants);
        final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());

        final ConsequenceOutput consequenceOutput = determiner.run(consensusPassedVariants);

        final List<SomaticVariant> potentialConsequentialMNVs =
                MNVDetector.locatePotentialMNVs(consensusPassedVariants, consequenceOutput.consequentialVariants());

        return new VariantAnalysis(variants, passedVariants, consensusPassedVariants, missenseVariants,
                consequenceOutput.consequentialVariants(), potentialConsequentialMNVs, consequenceOutput.findings());
    }

    @NotNull
    private static Predicate<SomaticVariant> isMissense() {
        return variant -> variant.hasConsequence(VariantConsequence.MISSENSE_VARIANT);
    }
}
