package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.isMissense;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.jetbrains.annotations.NotNull;

public class VariantInterpreter {

    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final ConsequenceRule consequenceRule;

    public static VariantInterpreter fromSlicingRegions(@NotNull final Slicer hmfSlicingRegion,
            @NotNull final Slicer giabHighConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        final ConsensusRule consensusRule = new ConsensusRule(giabHighConfidenceRegion, cpctSlicingRegion);
        final ConsequenceRule consequenceRule = new ConsequenceRule(hmfSlicingRegion, Lists.newArrayList());
        return new VariantInterpreter(consensusRule, consequenceRule);
    }

    private VariantInterpreter(@NotNull final ConsensusRule consensusRule,
            @NotNull final ConsequenceRule consequenceRule) {
        this.consensusRule = consensusRule;
        this.consequenceRule = consequenceRule;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final VCFSomaticFile variantFile) {
        final List<SomaticVariant> allVariants = variantFile.variants();
        final List<SomaticVariant> allPassedVariants = passOnly(allVariants);
        final List<SomaticVariant> consensusPassedVariants = consensusRule.apply(allPassedVariants);
        final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());
        final List<SomaticVariant> consequencePassedVariants = consequenceRule.apply(consensusPassedVariants);
        return new VariantAnalysis(allVariants, allPassedVariants, consensusPassedVariants, missenseVariants,
                consequencePassedVariants);
    }
}
