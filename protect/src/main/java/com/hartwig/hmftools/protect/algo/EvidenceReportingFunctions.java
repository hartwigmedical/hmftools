package com.hartwig.hmftools.protect.algo;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceComparator;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EvidenceReportingFunctions {

    private static final Set<Knowledgebase> TRIAL_SOURCES = Sets.newHashSet(Knowledgebase.ICLUSION, Knowledgebase.ACTIN);

    private EvidenceReportingFunctions() {
    }

    @NotNull
    public static List<ProtectEvidence> applyReportingAlgo(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> meetsMaxLevelSources = onlyReportWhenMeetsMaxLevelForSources(evidences);
        List<ProtectEvidence> maxLevelPerTreatmentEvent = onlyReportHighestLevelForTreatmentAndEvent(meetsMaxLevelSources);

        maxLevelPerTreatmentEvent.sort(new ProtectEvidenceComparator());
        return maxLevelPerTreatmentEvent;
    }

    @NotNull
    private static List<ProtectEvidence> onlyReportWhenMeetsMaxLevelForSources(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> result = Lists.newArrayList();

        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported()) {
                result.add(ImmutableProtectEvidence.builder()
                        .from(evidence)
                        .reported(meetsMaxReportableLevelForKnowledgebases(evidence))
                        .build());
            } else {
                result.add(evidence);
            }
        }
        return result;
    }

    private static boolean meetsMaxReportableLevelForKnowledgebases(@NotNull ProtectEvidence evidence) {
        EvidenceLevel lowestMaxReportingLevel = EvidenceLevel.A;
        for (ProtectSource source : evidence.sources()) {
            Knowledgebase knowledgebase = source.name();
            EvidenceLevel maxLevelForSource = evidence.direction().isCertain()
                    ? knowledgebase.maxCertainEvidenceReportingLevel()
                    : knowledgebase.maxPredictedEvidenceReportingLevel();

            if (lowestMaxReportingLevel.isHigher(maxLevelForSource)) {
                lowestMaxReportingLevel = maxLevelForSource;
            }
        }

        return !lowestMaxReportingLevel.isHigher(evidence.level());
    }

    @NotNull
    private static List<ProtectEvidence> onlyReportHighestLevelForTreatmentAndEvent(@NotNull List<ProtectEvidence> evidences) {
        Set<EvidenceKey> events = EvidenceKey.buildKeySet(evidences);

        List<ProtectEvidence> result = Lists.newArrayList();
        for (EvidenceKey event : events) {
            result.addAll(reportHighestPerEventTreatmentDirection(evidences.stream()
                    .filter(x -> x.direction().isResponsive())
                    .filter(x -> EvidenceKey.create(x).equals(event))
                    .collect(Collectors.toList())));

            result.addAll(reportHighestPerEventTreatmentDirection(evidences.stream()
                    .filter(x -> x.direction().isResistant())
                    .filter(x -> EvidenceKey.create(x).equals(event))
                    .collect(Collectors.toList())));

            result.addAll(reportHighestPerEventTreatmentDirection(evidences.stream()
                    .filter(x -> !x.direction().isResistant() &&  !x.direction().isResponsive())
                    .filter(x -> EvidenceKey.create(x).equals(event))
                    .collect(Collectors.toList())));
        }

        return result;
    }

    @NotNull
    private static List<ProtectEvidence> reportHighestPerEventTreatmentDirection(@NotNull Iterable<ProtectEvidence> evidences) {
        EvidenceLevel highestOnLabel = highestReportableLevel(true, evidences);
        EvidenceLevel highestOffLabel = highestReportableLevel(false, evidences);

        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            filtered.add(ImmutableProtectEvidence.builder()
                    .from(evidence)
                    .reported(reportEvidence(evidence, highestOnLabel, highestOffLabel))
                    .build());
        }

        return filtered;
    }

    private static boolean reportEvidence(@NotNull ProtectEvidence evidence, @Nullable EvidenceLevel highestOnLabel,
            @Nullable EvidenceLevel highestOffLabel) {
        if (evidence.reported()) {
            if (evidence.onLabel()) {
                assert highestOnLabel != null;
                return evidence.level() == highestOnLabel;
            } else {
                if (evidence.level() == highestOffLabel) {
                    return highestOnLabel == null || evidence.level().isHigher(highestOnLabel);
                }
            }
        }
        return false;
    }

    @Nullable
    @VisibleForTesting
    static EvidenceLevel highestReportableLevel(boolean isOnLabel, @NotNull Iterable<ProtectEvidence> evidences) {
        EvidenceLevel highest = null;
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported() && evidence.onLabel() == isOnLabel) {
                if (highest == null || evidence.level().isHigher(highest)) {
                    highest = evidence.level();
                }
            }
        }
        return highest;
    }

    @NotNull
    public static List<ProtectEvidence> reportOnLabelTrialsOnly(@NotNull Iterable<ProtectEvidence> evidences) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (isExclusiveTrialEvidence(evidence)) {
                result.add(ImmutableProtectEvidence.builder().from(evidence).reported(evidence.reported() && evidence.onLabel()).build());
            } else {
                result.add(evidence);
            }
        }
        return result;
    }

    private static boolean isExclusiveTrialEvidence(@NotNull ProtectEvidence evidence) {
        for (ProtectSource source : evidence.sources()) {
            if (!TRIAL_SOURCES.contains(source.name())) {
                return false;
            }
        }
        return true;
    }
}