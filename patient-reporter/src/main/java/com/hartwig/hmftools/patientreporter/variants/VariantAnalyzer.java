package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.isMissense;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.jetbrains.annotations.NotNull;

public class VariantAnalyzer {

    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final ConsequenceRule consequenceRule;

    public static VariantAnalyzer fromSlicingRegions(@NotNull final Slicer hmfSlicingRegion,
            @NotNull final Slicer giabHighConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        final ConsensusRule consensusRule = new ConsensusRule(giabHighConfidenceRegion, cpctSlicingRegion);
        final ConsequenceRule consequenceRule = ConsequenceRule.fromHmfSlicingRegion(hmfSlicingRegion);
        return new VariantAnalyzer(consensusRule, consequenceRule);
    }

    private VariantAnalyzer(@NotNull final ConsensusRule consensusRule,
            @NotNull final ConsequenceRule consequenceRule) {
        this.consensusRule = consensusRule;
        this.consequenceRule = consequenceRule;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> passedVariants = passOnly(variants);
        final List<SomaticVariant> consensusPassedVariants = consensusRule.apply(passedVariants);
        final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());
        final List<SomaticVariant> consequencePassedVariants = consequenceRule.apply(consensusPassedVariants);
        return new VariantAnalysis(passedVariants, consensusPassedVariants, missenseVariants,
                consequencePassedVariants);
    }
}
