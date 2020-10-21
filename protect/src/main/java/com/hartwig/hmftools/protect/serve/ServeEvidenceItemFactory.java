package com.hartwig.hmftools.protect.serve;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public class ServeEvidenceItemFactory {

    @NotNull
    public static ServeEvidenceItem create(@NotNull String genomicEvent, @NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return builder(doid, actionable).genomicEvent(genomicEvent).reported(true).build();
    }

    @NotNull
    public static ImmutableServeEvidenceItem.Builder builder(@NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return ImmutableServeEvidenceItem.builder()
                .source(actionable.source())
                .direction(actionable.direction())
                .level(actionable.level())
                .treatment(actionable.treatment())
                .onLabel(doid.contains(actionable.doid()))
                .direction(actionable.direction())
                .url(actionable.url());
    }

    @NotNull
    public static List<ServeEvidenceItem> doNotReportInsignificantEvidence(@NotNull Collection<ServeEvidenceItem> evidence) {
        final EvidenceLevel highestResponse = minLevel(EvidenceDirection.RESPONSIVE, evidence);
        final EvidenceLevel highestResistance = minLevel(EvidenceDirection.RESISTANT, evidence);

        return evidence.stream()
                .map(x -> ImmutableServeEvidenceItem.builder()
                        .from(x)
                        .reported(x.reported() && report(x, highestResponse, highestResistance))
                        .build())
                .collect(Collectors.toList());
    }

    static boolean report(ServeEvidenceItem victim, EvidenceLevel minResponse, EvidenceLevel minResistance) {
        return victim.direction().equals(EvidenceDirection.RESPONSIVE)
                ? victim.level().ordinal() <= minResponse.ordinal()
                : victim.level().ordinal() <= minResistance.ordinal();
    }

    @NotNull
    static EvidenceLevel minLevel(@NotNull EvidenceDirection direction, @NotNull Collection<? extends ServeEvidenceItem> actionable) {
        Optional<EvidenceLevel> highestOnLabel = actionable.stream()
                .filter(ServeEvidenceItem::onLabel)
                .filter(x -> x.direction().equals(direction))
                .min(Comparator.comparing(ServeEvidenceItem::level))
                .map(ServeEvidenceItem::level);

        return highestOnLabel.orElseGet(() -> actionable.stream()
                .filter(x -> !x.onLabel())
                .filter(x -> x.direction().equals(direction))
                .min(Comparator.comparing(ServeEvidenceItem::level))
                .map(ServeEvidenceItem::level)
                .orElse(EvidenceLevel.D));
    }

}
