package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.protect.variants.DriverInterpretation;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;
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
    private final List<ActionableHotspot> hotspots;
    @NotNull
    private final List<ActionableRange> ranges;
    @NotNull
    private final List<ActionableGene> genes;

    public VariantEvidence(@NotNull final List<ActionableHotspot> hotspots, @NotNull final List<ActionableRange> ranges,
            @NotNull final List<ActionableGene> genes) {
        this.hotspots = hotspots;
        this.ranges = ranges;
        this.genes = genes.stream()
                .filter(x -> x.event() == GeneLevelEvent.ACTIVATION || x.event() == GeneLevelEvent.INACTIVATION
                        || x.event() == GeneLevelEvent.ANY_MUTATION)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull List<ReportableVariant> germline,
            @NotNull List<ReportableVariant> somatic) {
        List<ReportableVariant> variants = ReportableVariantFactory.mergeVariantLists(germline, somatic);
        return variants.stream().flatMap(x -> evidence(doids, x).stream()).collect(Collectors.toList());
    }

    @NotNull
    private List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull ReportableVariant reportable) {
        List<ProtectEvidenceItem> hotspotEvidence = hotspots.stream()
                .filter(x -> hotspotMatch(x, reportable))
                .map(x -> evidence(true, doids, reportable, x))
                .collect(Collectors.toList());

        List<ProtectEvidenceItem> rangeEvidence = ranges.stream()
                .filter(x -> rangeMatch(x, reportable))
                .map(x -> evidence(true, doids, reportable, x))
                .collect(Collectors.toList());

        List<ProtectEvidenceItem> geneEvidence = genes.stream()
                .filter(x -> geneMatch(x, reportable))
                .map(x -> evidence(reportable.driverLikelihoodInterpretation() == DriverInterpretation.HIGH, doids, reportable, x))
                .collect(Collectors.toList());

        List<ProtectEvidenceItem> result = Lists.newArrayList();
        result.addAll(hotspotEvidence);
        result.addAll(rangeEvidence);
        result.addAll(geneEvidence);

        return ProtectEvidenceItems.reportHighest(result);
    }

    private static boolean hotspotMatch(@NotNull ActionableHotspot hotspot, @NotNull ReportableVariant variant) {
        return variant.chromosome().equals(hotspot.chromosome()) && hotspot.position() == variant.position() && hotspot.ref()
                .equals(variant.ref()) && hotspot.alt().equals(hotspot.alt());
    }

    private static boolean rangeMatch(@NotNull ActionableRange range, @NotNull ReportableVariant variant) {
        return variant.chromosome().equals(range.chromosome()) && variant.gene().equals(range.gene()) && variant.position() >= range.start()
                && variant.position() <= range.end() && meetsMutationTypeFilter(range.mutationType(), variant);
    }

    private static boolean meetsMutationTypeFilter(@NotNull MutationTypeFilter filter, @NotNull ReportableVariant variant) {
        CodingEffect effect = variant.canonicalCodingEffect();
        switch (filter) {
            case NONSENSE_OR_FRAMESHIFT:
                return effect == CodingEffect.NONSENSE_OR_FRAMESHIFT;
            case SPLICE:
                return effect == CodingEffect.SPLICE;
            case INFRAME:
                return effect == CodingEffect.MISSENSE && (isInsert(variant) || isDelete(variant));
            case INFRAME_DELETION:
                return effect == CodingEffect.MISSENSE && isDelete(variant);
            case INFRAME_INSERTION:
                return effect == CodingEffect.MISSENSE && isInsert(variant);
            case MISSENSE:
                return effect == CodingEffect.MISSENSE;
            default: {
                LOGGER.warn("Unrecognized mutation type filter: '{}'", filter);
                return false;
            }
        }
    }

    private static boolean isInsert(@NotNull ReportableVariant variant) {
        return variant.alt().length() > variant.ref().length();
    }

    private static boolean isDelete(@NotNull ReportableVariant variant) {
        return variant.alt().length() < variant.ref().length();
    }

    private static boolean geneMatch(@NotNull ActionableGene gene, @NotNull ReportableVariant variant) {
        assert gene.event() == GeneLevelEvent.ACTIVATION || gene.event() == GeneLevelEvent.INACTIVATION
                || gene.event() == GeneLevelEvent.ANY_MUTATION;

        return gene.gene().equals(variant.gene());
    }

    @NotNull
    private static ProtectEvidenceItem evidence(boolean report, @NotNull Set<String> doids, @NotNull ReportableVariant reportable,
            @NotNull ActionableEvent actionable) {
        return ProtectEvidenceItems.builder(doids, actionable)
                .genomicEvent(reportable.genomicEvent())
                .germline(reportable.source() == ReportableVariantSource.GERMLINE)
                .reported(report)
                .build();
    }
}
