package com.hartwig.hmftools.protect.variants;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.protect.serve.ServeEvidenceItem;
import com.hartwig.hmftools.protect.serve.ServeEvidenceItemFactory;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;

import org.jetbrains.annotations.NotNull;

public class ReportableVariantEvidence {

    private static final Set<CodingEffect> RANGE_CODING_EFFECTS =
            Sets.newHashSet(CodingEffect.SPLICE, CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE);

    private final List<ActionableRange> ranges;
    private final List<ActionableHotspot> hotspots;

    public ReportableVariantEvidence(@NotNull final List<ActionableHotspot> hotspots, @NotNull final List<ActionableRange> ranges) {
        this.hotspots = hotspots;
        this.ranges = ranges;
    }

    @NotNull
    public List<ServeEvidenceItem> evidence(@NotNull Set<String> doid, @NotNull List<ReportableVariant> germline,
            @NotNull List<ReportableVariant> somatic) {
        List<ReportableVariant> variants = ReportableVariantFactory.mergeSomaticAndGermlineVariants(germline, somatic);
        return variants.stream().flatMap(x -> evidence(doid, x).stream()).collect(Collectors.toList());
    }

    @NotNull
    public List<ServeEvidenceItem> evidence(@NotNull Set<String> doid, @NotNull ReportableVariant reportable) {
        boolean report = reportable.driverLikelihoodInterpretation().equals(DriverInterpretation.HIGH);

        List<ServeEvidenceItem> hotspotEvidence = hotspots.stream()
                .filter(x -> hotspotMatch(x, reportable))
                .map(x -> evidence(report, doid, reportable, x))
                .collect(Collectors.toList());

        List<ServeEvidenceItem> rangeEvidence = ranges.stream()
                .filter(x -> rangeMatch(x, reportable))
                .map(x -> evidence(report, doid, reportable, x))
                .collect(Collectors.toList());

        Set<ServeEvidenceItem> result = Sets.newHashSet();
        result.addAll(hotspotEvidence);
        result.addAll(rangeEvidence);

        return ServeEvidenceItemFactory.doNotReportInsignificantEvidence(result);
    }

    private static boolean rangeMatch(ActionableRange range, ReportableVariant variant) {
        return RANGE_CODING_EFFECTS.contains(variant.canonicalCodingEffect()) && variant.chromosome().equals(range.chromosome())
                && range.gene().equals(range.gene()) && variant.position() >= range.start() && variant.position() <= range.end();
    }

    private static boolean hotspotMatch(ActionableHotspot hotspot, ReportableVariant variant) {
        return variant.chromosome().equals(hotspot.chromosome()) && hotspot.alt().equals(hotspot.alt())
                && hotspot.position() == variant.position() && hotspot.ref().equals(variant.ref());
    }

    @NotNull
    private static ServeEvidenceItem evidence(boolean report, @NotNull Set<String> doid, @NotNull ReportableVariant reportable,
            @NotNull ActionableEvent actionable) {
        return ServeEvidenceItemFactory.builder(doid, actionable).genomicEvent(reportable.genomicEvent()).reported(report).build();
    }

}
