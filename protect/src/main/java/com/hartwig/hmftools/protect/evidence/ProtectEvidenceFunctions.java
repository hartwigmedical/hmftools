package com.hartwig.hmftools.protect.evidence;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceFunctions {

    private ProtectEvidenceFunctions() {
    }

    @NotNull
    public static ImmutableProtectEvidence.Builder builder(@NotNull Set<String> doids, @NotNull ActionableEvent actionable) {
        return ImmutableProtectEvidence.builder()
                .source(actionable.source())
                .treatment(actionable.treatment())
                .level(actionable.level())
                .direction(actionable.direction())
                .onLabel(doids.contains(actionable.doid()))
                .urls(actionable.urls());
    }

    @NotNull
    public static List<ProtectEvidence> reportHighest(@NotNull Collection<ProtectEvidence> evidence) {
        final Set<String> events = evidence.stream().map(ProtectEvidence::genomicEvent).collect(Collectors.toSet());
        final Set<String> treatments = evidence.stream().map(ProtectEvidence::treatment).collect(Collectors.toSet());

        List<ProtectEvidence> result = Lists.newArrayList();
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
    private static List<ProtectEvidence> reportHighestPerEventTreatmentDirection(@NotNull Collection<ProtectEvidence> evidence) {
        Optional<EvidenceLevel> highestOnLabel = highestReportableLevel(true, evidence);
        Optional<EvidenceLevel> highestOffLabel = highestReportableLevel(false, evidence);
        Predicate<ProtectEvidence> report = x -> x.reported() && x.onLabel()
                ? highestOnLabel.filter(highest -> x.level().ordinal() == highest.ordinal()).isPresent()
                : highestOffLabel.filter(highest -> x.level().ordinal() == highest.ordinal()).isPresent()
                        && !highestOnLabel.filter(highest -> x.level().ordinal() >= highest.ordinal()).isPresent();

        return evidence.stream()
                .map(x -> ImmutableProtectEvidence.builder().from(x).reported(report.test(x)).build())
                .collect(Collectors.toList());
    }

    @NotNull
    @VisibleForTesting
    static Optional<EvidenceLevel> highestReportableLevel(boolean isOnLabel,
            @NotNull Collection<? extends ProtectEvidence> actionable) {
        return actionable.stream()
                .filter(x -> x.level().ordinal() <= EvidenceLevel.B.ordinal())
                .filter(x -> x.onLabel() == isOnLabel)
                .filter(x -> x.reported())
                .min(Comparator.comparing(ProtectEvidence::level))
                .map(ProtectEvidence::level);
    }
}
