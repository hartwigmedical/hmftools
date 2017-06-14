package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.HmfSlicer;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;

import org.jetbrains.annotations.NotNull;

public class VariantAnalyzer {

    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final ConsequenceDeterminer determiner;

    public static VariantAnalyzer fromSlicingRegions(@NotNull final HmfSlicer hmfSlicingRegion,
            @NotNull final Slicer giabHighConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        final ConsensusRule consensusRule = ConsensusRule.fromSlicers(giabHighConfidenceRegion, cpctSlicingRegion);
        final ConsequenceDeterminer determiner = fromHmfSlicingRegion(hmfSlicingRegion);
        return new VariantAnalyzer(consensusRule, determiner);
    }

    @NotNull
    private static ConsequenceDeterminer fromHmfSlicingRegion(@NotNull final HmfSlicer hmfSlicingRegion) {
        return new ConsequenceDeterminer(hmfSlicingRegion, extractTranscriptMap(hmfSlicingRegion));
    }

    @NotNull
    private static Map<String, HmfGenomeRegion> extractTranscriptMap(final @NotNull HmfSlicer hmfSlicingRegion) {
        final Map<String, HmfGenomeRegion> transcriptMap = Maps.newHashMap();
        for (final HmfGenomeRegion region : hmfSlicingRegion.hmfRegions()) {
            transcriptMap.put(region.transcriptID(), region);
        }
        return transcriptMap;
    }

    private VariantAnalyzer(@NotNull final ConsensusRule consensusRule,
            @NotNull final ConsequenceDeterminer determiner) {
        this.consensusRule = consensusRule;
        this.determiner = determiner;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> passedVariants = passOnly(variants);
        final List<SomaticVariant> consensusPassedVariants = consensusRule.removeUnreliableVariants(passedVariants);
        final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());

        final ConsequenceOutput consequenceOutput = determiner.run(consensusPassedVariants);

        final List<SomaticVariant> potentialConsequentialMNVs = MNVDetector.locatePotentialMNVs(
                consensusPassedVariants, consequenceOutput.consequentialVariants());

        return new VariantAnalysis(variants, passedVariants, consensusPassedVariants, missenseVariants,
                consequenceOutput.consequentialVariants(), potentialConsequentialMNVs, consequenceOutput.findings());
    }

    @NotNull
    private static Predicate<SomaticVariant> isMissense() {
        return variant -> variant.hasConsequence(VariantConsequence.MISSENSE_VARIANT);
    }
}
