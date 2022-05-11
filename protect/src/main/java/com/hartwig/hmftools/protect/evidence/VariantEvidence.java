package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VariantEvidence {

    private static final Logger LOGGER = LogManager.getLogger(VariantEvidence.class);

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableHotspot> hotspots;
    @NotNull
    private final List<ActionableRange> ranges;
    @NotNull
    private final List<ActionableGene> genes;

    public VariantEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableHotspot> hotspots, @NotNull final List<ActionableRange> ranges,
            @NotNull final List<ActionableGene> genes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.hotspots = hotspots;
        this.ranges = ranges;
        this.genes = genes.stream()
                .filter(x -> x.event() == GeneLevelEvent.ACTIVATION || x.event() == GeneLevelEvent.INACTIVATION
                        || x.event() == GeneLevelEvent.ANY_MUTATION)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<ReportableVariant> germline, @NotNull List<ReportableVariant> somatic,
            @NotNull List<SomaticVariant> unreportedSomaticVariants) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (ReportableVariant reportableVariant : ReportableVariantFactory.mergeVariantLists(germline, somatic)) {
            evidences.addAll(evidence(reportableVariant,
                    reportableVariant.driverLikelihoodInterpretation(),
                    reportableVariant.source() == ReportableVariantSource.GERMLINE,
                    true,
                    reportableVariant.transcript(),
                    reportableVariant.isCanonical(),
                    reportableVariant.otherReportedEffects()));
        }

        for (SomaticVariant unreportedVariant : unreportedSomaticVariants) {
            evidences.addAll(evidence(unreportedVariant,
                    DriverInterpretation.LOW,
                    false,
                    false,
                    unreportedVariant.canonicalTranscript(),
                    true,
                    unreportedVariant.otherReportedEffects()));
        }

        return evidences;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull Variant variant, @NotNull DriverInterpretation driverInterpretation, boolean germline,
            boolean mayReport, @NotNull String transcript, boolean isCanonical, @NotNull String otherReportedEffects) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (ActionableHotspot hotspot : hotspots) {
            if (hotspotMatch(variant, hotspot)) {
                evidences.add(evidence(variant,
                        germline,
                        mayReport,
                        hotspot,
                        driverInterpretation,
                        transcript,
                        isCanonical,
                        otherReportedEffects));
            }
        }

        for (ActionableRange range : ranges) {
            if (rangeMatch(variant, range, isCanonical, otherReportedEffects)) {
                evidences.add(evidence(variant,
                        germline,
                        mayReport,
                        range,
                        driverInterpretation,
                        transcript,
                        isCanonical,
                        otherReportedEffects));
            }
        }

        for (ActionableGene gene : genes) {
            if (geneMatch(variant, gene, isCanonical, otherReportedEffects)) {
                evidences.add(evidence(variant,
                        germline,
                        mayReport && driverInterpretation == DriverInterpretation.HIGH,
                        gene,
                        driverInterpretation,
                        transcript,
                        isCanonical,
                        otherReportedEffects));
            }
        }

        return evidences;
    }

    @NotNull
    private ProtectEvidence evidence(@NotNull Variant variant, boolean germline, boolean report, @NotNull ActionableEvent actionable,
            @NotNull DriverInterpretation driverInterpretation, @NotNull String transcript, boolean isCanonical,
            @NotNull String otherReportedEffects) {

        return personalizedEvidenceFactory.evidenceBuilder(actionable)
                .gene(variant.gene())
                .transcript(transcript)
                .isCanonical(isCanonical)
                .event(determineEvent(isCanonical, variant, otherReportedEffects))
                .germline(germline)
                .reported(report)
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretVariant(driverInterpretation))
                .build();
    }

    @VisibleForTesting
    static String determineEvent(boolean isCanonical, @NotNull Variant variant, @NotNull String otherReportedEffects) {
        return isCanonical
                ? ProtectEventGenerator.variantEvent(variant)
                : ProtectEventGenerator.variantEventNonCanonical(otherReportedEffects);
    }

    private static boolean hotspotMatch(@NotNull Variant variant, @NotNull ActionableHotspot hotspot) {
        return variant.chromosome().equals(hotspot.chromosome()) && hotspot.position() == variant.position() && hotspot.ref()
                .equals(variant.ref()) && hotspot.alt().equals(variant.alt());
    }

    private static boolean rangeMatch(@NotNull Variant variant, @NotNull ActionableRange range, boolean isCanonical,
            @NotNull String otherReportedEffects) {
        return variant.chromosome().equals(range.chromosome()) && variant.gene().equals(range.gene()) && variant.position() >= range.start()
                && variant.position() <= range.end() && meetsMutationTypeFilter(variant,
                range.mutationType(),
                isCanonical,
                otherReportedEffects);
    }

    private static boolean geneMatch(@NotNull Variant variant, @NotNull ActionableGene gene, boolean isCanonical,
            @NotNull String otherReportedEffects) {
        assert gene.event() == GeneLevelEvent.ACTIVATION || gene.event() == GeneLevelEvent.INACTIVATION
                || gene.event() == GeneLevelEvent.ANY_MUTATION;

        return gene.gene().equals(variant.gene()) && meetsMutationTypeFilter(variant,
                MutationTypeFilter.ANY,
                isCanonical,
                otherReportedEffects);
    }

    @NotNull
    @VisibleForTesting
    static CodingEffect extractCodingEffectOther(@NotNull String otherReportedEffects) {
        return !otherReportedEffects.isEmpty() ? CodingEffect.valueOf(otherReportedEffects.split("\\|")[4]) : CodingEffect.UNDEFINED;
    }

    private static boolean meetsMutationTypeFilter(@NotNull Variant variant, @NotNull MutationTypeFilter filter, boolean isCanonical,
            @NotNull String otherReportedEffects) {
        CodingEffect effectOther = extractCodingEffectOther(otherReportedEffects);
        CodingEffect effectCanonical = variant.canonicalCodingEffect();
        CodingEffect effect = isCanonical ? effectCanonical : effectOther;

        switch (filter) {
            case NONSENSE_OR_FRAMESHIFT:
                return effect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
            case SPLICE:
                return effect == CodingEffect.SPLICE;
            case INFRAME:
                return effect == CodingEffect.MISSENSE && variant.type() == VariantType.INDEL;
            case INFRAME_DELETION:
                return effect == CodingEffect.MISSENSE && isDelete(variant);
            case INFRAME_INSERTION:
                return effect == CodingEffect.MISSENSE && isInsert(variant);
            case MISSENSE:
                return effect == CodingEffect.MISSENSE;
            case ANY:
                return effect == CodingEffect.MISSENSE || effect == CodingEffect.NONSENSE_OR_FRAMESHIFT || effect == CodingEffect.SPLICE;
            default: {
                LOGGER.warn("Unrecognized mutation type filter: '{}'", filter);
                return false;
            }
        }
    }

    private static boolean isInsert(@NotNull Variant variant) {
        return variant.type() == VariantType.INDEL && variant.alt().length() > variant.ref().length();
    }

    private static boolean isDelete(@NotNull Variant variant) {
        return variant.type() == VariantType.INDEL && variant.alt().length() < variant.ref().length();
    }
}
