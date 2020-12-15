package com.hartwig.hmftools.protect.evidence;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceItems {

    private ProtectEvidenceItems() {
    }

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
                .urls(actionable.urls());
    }

    @NotNull
    public static List<ProtectEvidenceItem> reportHighest(@NotNull Collection<ProtectEvidenceItem> evidence) {
        final Set<String> events = evidence.stream().map(ProtectEvidenceItem::genomicEvent).collect(Collectors.toSet());
        final Set<String> treatments = evidence.stream().map(ProtectEvidenceItem::treatment).collect(Collectors.toSet());

        List<ProtectEvidenceItem> result = Lists.newArrayList();
        for (String event : events) {
            for (String treatment : treatments) {
                for (EvidenceDirection direction : EvidenceDirection.values()) {
                    result.addAll(reportHighestPerEventTreatmentDirection(evidence.stream()
                            .filter(x -> x.treatment().equals(treatment))
                            .filter(x -> x.direction().equals(direction))
                            .filter(x -> x.genomicEvent().equals(event))
                            .collect(Collectors.toSet())));
                }
            }
        }

        return result.stream().sorted().collect(Collectors.toList());
    }

    @NotNull
    private static List<ProtectEvidenceItem> reportHighestPerEventTreatmentDirection(@NotNull Collection<ProtectEvidenceItem> evidence) {
        Optional<EvidenceLevel> highestOnLabel = highestReportableLevel(true, evidence);
        Optional<EvidenceLevel> highestOffLabel = highestReportableLevel(false, evidence);
        Predicate<ProtectEvidenceItem> report = x -> x.reported() && x.onLabel()
                ? highestOnLabel.filter(highest -> x.level().ordinal() == highest.ordinal()).isPresent()
                : highestOffLabel.filter(highest -> x.level().ordinal() == highest.ordinal()).isPresent()
                        && !highestOnLabel.filter(highest -> x.level().ordinal() >= highest.ordinal()).isPresent();

        return evidence.stream()
                .map(x -> ImmutableProtectEvidenceItem.builder().from(x).reported(report.test(x)).build())
                .collect(Collectors.toList());

    }

    @NotNull
    static Optional<EvidenceLevel> highestReportableLevel(boolean isOnLabel,
            @NotNull Collection<? extends ProtectEvidenceItem> actionable) {
        return actionable.stream()
                .filter(x -> x.level().ordinal() <= EvidenceLevel.B.ordinal())
                .filter(x -> x.onLabel() == isOnLabel)
                .min(Comparator.comparing(ProtectEvidenceItem::level))
                .map(ProtectEvidenceItem::level);
    }
}
