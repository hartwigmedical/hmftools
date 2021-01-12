package com.hartwig.hmftools.protect.evidence;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public final class ProtectEvidenceFunctions {

    private ProtectEvidenceFunctions() {
    }

    @NotNull
    public static List<ProtectEvidence> consolidate(@NotNull List<ProtectEvidence> evidences) {
        Map<ProtectEvidence, ConsolidatedData> dataPerEvidence = Maps.newHashMap();
        for (ProtectEvidence evidence : evidences) {
            ProtectEvidence strippedEvidence =
                    ImmutableProtectEvidence.builder().from(evidence).urls(Sets.newHashSet()).sources(Sets.newHashSet()).build();
            ConsolidatedData data = dataPerEvidence.get(strippedEvidence);
            if (data == null) {
                data = new ConsolidatedData();
            }
            data.appendEvidence(evidence);
            dataPerEvidence.put(strippedEvidence, data);
        }

        List<ProtectEvidence> consolidatedEvents = Lists.newArrayList();
        for (Map.Entry<ProtectEvidence, ConsolidatedData> entry : dataPerEvidence.entrySet()) {
            consolidatedEvents.add(ImmutableProtectEvidence.builder()
                    .from(entry.getKey())
                    .sources(entry.getValue().sources())
                    .urls(entry.getValue().urls())
                    .build());
        }
        return consolidatedEvents;
    }

    @NotNull
    public static List<ProtectEvidence> reportHighest(@NotNull List<ProtectEvidence> evidence) {
        Set<String> events = evidence.stream().map(ProtectEvidence::genomicEvent).collect(Collectors.toSet());
        Set<String> treatments = evidence.stream().map(ProtectEvidence::treatment).collect(Collectors.toSet());

        List<ProtectEvidence> result = Lists.newArrayList();
        for (String event : events) {
            for (String treatment : treatments) {
                for (EvidenceDirection direction : EvidenceDirection.values()) {
                    result.addAll(reportHighestPerEventTreatmentDirection(evidence.stream()
                            .filter(x -> x.treatment().equals(treatment))
                            .filter(x -> x.direction().equals(direction))
                            .filter(x -> x.genomicEvent().equals(event))
                            .collect(Collectors.toList())));
                }
            }
        }

        return result.stream().sorted().collect(Collectors.toList());
    }

    @NotNull
    private static List<ProtectEvidence> reportHighestPerEventTreatmentDirection(@NotNull List<ProtectEvidence> evidence) {
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
    static Optional<EvidenceLevel> highestReportableLevel(boolean isOnLabel, @NotNull List<ProtectEvidence> actionable) {
        return actionable.stream()
                .filter(x -> x.level().ordinal() <= EvidenceLevel.B.ordinal())
                .filter(x -> x.onLabel() == isOnLabel)
                .filter(x -> x.reported())
                .min(Comparator.comparing(ProtectEvidence::level))
                .map(ProtectEvidence::level);
    }

    private static class ConsolidatedData {

        @NotNull
        private final Set<Knowledgebase> sources = Sets.newHashSet();
        @NotNull
        private final Set<String> urls = Sets.newHashSet();

        public ConsolidatedData() {
        }

        public void appendEvidence(@NotNull ProtectEvidence evidence) {
            sources.addAll(evidence.sources());
            urls.addAll(evidence.urls());
        }

        @NotNull
        public Set<Knowledgebase> sources() {
            return sources;
        }

        @NotNull
        public Set<String> urls() {
            return urls;
        }
    }
}
