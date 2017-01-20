package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(VariantAnalyzer.class);

    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final Reporter reporter;

    public static VariantAnalyzer fromSlicingRegions(@NotNull final Slicer hmfSlicingRegion,
            @NotNull final Slicer giabHighConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        final ConsensusRule consensusRule = new ConsensusRule(giabHighConfidenceRegion, cpctSlicingRegion);
        final Reporter reporter = fromHmfSlicingRegion(hmfSlicingRegion);
        return new VariantAnalyzer(consensusRule, reporter);
    }

    @NotNull
    private static Reporter fromHmfSlicingRegion(@NotNull final Slicer hmfSlicingRegion) {
        return new Reporter(hmfSlicingRegion, extractTranscriptMap(hmfSlicingRegion));
    }

    @NotNull
    private static Map<String, HMFSlicingAnnotation> extractTranscriptMap(final @NotNull Slicer hmfSlicingRegion) {
        final Map<String, HMFSlicingAnnotation> transcriptMap = Maps.newHashMap();
        for (final GenomeRegion region : hmfSlicingRegion.regions()) {
            final HMFSlicingAnnotation annotation = HMFSlicingAnnotation.fromGenomeRegion(region);
            if (annotation == null) {
                LOGGER.warn("Could not extract annotation from hmf slicing region: " + region);
            } else {
                transcriptMap.put(annotation.transcriptID(), annotation);
            }
        }
        return transcriptMap;
    }

    private VariantAnalyzer(@NotNull final ConsensusRule consensusRule, @NotNull final Reporter reporter) {
        this.consensusRule = consensusRule;
        this.reporter = reporter;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> passedVariants = passOnly(variants);
        final List<SomaticVariant> consensusPassedVariants = consensusRule.apply(passedVariants);
        final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());

        final List<VariantReport> variantsToReport = reporter.run(consensusPassedVariants);

        return new VariantAnalysis(passedVariants, consensusPassedVariants, missenseVariants, variantsToReport);
    }

    @NotNull
    private static Predicate<SomaticVariant> isMissense() {
        return variant -> variant.hasConsequence(VariantConsequence.MISSENSE_VARIANT);
    }
}
