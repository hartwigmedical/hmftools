package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.protect.variants.DriverInterpretation;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantFactory;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;

import org.jetbrains.annotations.NotNull;

public class VariantEvidence {

    private static final Set<CodingEffect> RANGE_CODING_EFFECTS =
            Sets.newHashSet(CodingEffect.SPLICE, CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE);

    @NotNull
    private final List<ActionableHotspot> hotspots;
    @NotNull
    private final List<ActionableRange> ranges;

    public VariantEvidence(@NotNull final List<ActionableHotspot> hotspots, @NotNull final List<ActionableRange> ranges) {
        this.hotspots = hotspots;
        this.ranges = ranges;
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull List<ReportableVariant> germline,
            @NotNull List<ReportableVariant> somatic) {
        List<ReportableVariant> variants = ReportableVariantFactory.mergeVariantLists(germline, somatic);
        return variants.stream().flatMap(x -> evidence(doids, x).stream()).collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull ReportableVariant reportable) {
        boolean report = reportable.driverLikelihoodInterpretation().equals(DriverInterpretation.HIGH);

        List<ProtectEvidenceItem> hotspotEvidence = hotspots.stream()
                .filter(x -> hotspotMatch(x, reportable))
                .map(x -> evidence(true, doids, reportable, x))
                .collect(Collectors.toList());

        // TODO (DEV-642) Match mutation type against actionable mutation type filter
        List<ProtectEvidenceItem> rangeEvidence = ranges.stream()
                .filter(x -> rangeMatch(x, reportable))
                .map(x -> evidence(report, doids, reportable, x))
                .collect(Collectors.toList());

        // TODO (DEV-642) Include match for INACTIVATION/ACTIVATION/ANY_MUTATION

        Set<ProtectEvidenceItem> result = Sets.newHashSet();
        result.addAll(hotspotEvidence);
        result.addAll(rangeEvidence);

        return ProtectEvidenceItems.reportHighest(result);
    }

    private static boolean rangeMatch(@NotNull ActionableRange range, @NotNull ReportableVariant variant) {
        return RANGE_CODING_EFFECTS.contains(variant.canonicalCodingEffect()) && variant.chromosome().equals(range.chromosome())
                && range.gene().equals(range.gene()) && variant.position() >= range.start() && variant.position() <= range.end();
    }

    private static boolean hotspotMatch(@NotNull ActionableHotspot hotspot, @NotNull ReportableVariant variant) {
        return variant.chromosome().equals(hotspot.chromosome()) && hotspot.alt().equals(hotspot.alt())
                && hotspot.position() == variant.position() && hotspot.ref().equals(variant.ref());
    }

    @NotNull
    private static ProtectEvidenceItem evidence(boolean report, @NotNull Set<String> doids, @NotNull ReportableVariant reportable,
            @NotNull ActionableEvent actionable) {
        return ProtectEvidenceItems.builder(doids, actionable).genomicEvent(reportable.genomicEvent()).reported(report).build();
    }
}
