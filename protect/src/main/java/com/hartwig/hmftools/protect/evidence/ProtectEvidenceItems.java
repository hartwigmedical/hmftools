package com.hartwig.hmftools.protect.evidence;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.protect.ImmutableProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.serve.EvidenceDirection;
import com.hartwig.hmftools.common.serve.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public class ProtectEvidenceItems {

    @NotNull
    public static ProtectEvidenceItem create(@NotNull String genomicEvent, @NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return builder(doid, actionable).genomicEvent(genomicEvent).reported(true).build();
    }

    @NotNull
    public static ImmutableProtectEvidenceItem.Builder builder(@NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return ImmutableProtectEvidenceItem.builder()
                .source(actionable.source())
                .direction(actionable.direction())
                .level(actionable.level())
                .treatment(actionable.treatment())
                .onLabel(doid.contains(actionable.doid()))
                .direction(actionable.direction())
                .url(actionable.url());
    }

    @NotNull
    public static List<ProtectEvidenceItem> doNotReportInsignificantEvidence(@NotNull Collection<ProtectEvidenceItem> evidence) {
        final EvidenceLevel highestResponse = highestLevel(EvidenceDirection.RESPONSIVE, evidence);
        final EvidenceLevel highestResistance = highestLevel(EvidenceDirection.RESISTANT, evidence);

        return evidence.stream()
                .map(x -> ImmutableProtectEvidenceItem.builder()
                        .from(x)
                        .reported(x.reported() && report(x, highestResponse, highestResistance))
                        .build())
                .collect(Collectors.toList());
    }

    static boolean report(ProtectEvidenceItem victim, EvidenceLevel minResponse, EvidenceLevel minResistance) {
        return victim.direction().equals(EvidenceDirection.RESPONSIVE)
                ? victim.level().ordinal() <= minResponse.ordinal()
                : victim.level().ordinal() <= minResistance.ordinal();
    }

    @NotNull
    static EvidenceLevel highestLevel(@NotNull EvidenceDirection direction, @NotNull Collection<? extends ProtectEvidenceItem> actionable) {
        Optional<EvidenceLevel> highestOnLabel = actionable.stream()
                .filter(ProtectEvidenceItem::onLabel)
                .filter(x -> x.direction().equals(direction))
                .min(Comparator.comparing(ProtectEvidenceItem::level))
                .map(ProtectEvidenceItem::level);

        return highestOnLabel.orElseGet(() -> actionable.stream()
                .filter(x -> !x.onLabel())
                .filter(x -> x.direction().equals(direction))
                .min(Comparator.comparing(ProtectEvidenceItem::level))
                .map(ProtectEvidenceItem::level)
                .orElse(EvidenceLevel.D));
    }

}
