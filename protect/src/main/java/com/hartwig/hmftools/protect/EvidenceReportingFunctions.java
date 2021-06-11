package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EvidenceReportingFunctions {

    private static final Set<Knowledgebase> TRIAL_SOURCES = Sets.newHashSet(Knowledgebase.ICLUSION);

    private EvidenceReportingFunctions() {
    }

    @NotNull
    public static List<ProtectEvidence> applyReportingAlgo(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> meetsMaxSources = onlyReportWhenMeetsMaxLevelForSources(evidences);
        List<ProtectEvidence> maxLevelPerTreatmentEvent = onlyReportHighestLevelForTreatmentAndEvent(meetsMaxSources);

        return maxLevelPerTreatmentEvent.stream().sorted().collect(Collectors.toList());
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
        for (Knowledgebase source : evidence.sources()) {
            EvidenceLevel maxLevelForSource = evidence.direction().isCertain()
                    ? source.maxCertainEvidenceReportingLevel()
                    : source.maxCertainEvidenceReportingLevel();

            if (lowestMaxReportingLevel.isHigher(maxLevelForSource)) {
                lowestMaxReportingLevel = source.maxCertainEvidenceReportingLevel();
            }
        }

        return !lowestMaxReportingLevel.isHigher(evidence.level());
    }

    @NotNull
    private static List<ProtectEvidence> onlyReportHighestLevelForTreatmentAndEvent(@NotNull List<ProtectEvidence> evidences) {
        Set<String> events = evidences.stream().map(ProtectEvidence::genomicEvent).collect(Collectors.toSet());
        Set<String> treatments = evidences.stream().map(ProtectEvidence::treatment).collect(Collectors.toSet());

        List<ProtectEvidence> result = Lists.newArrayList();
        for (String event : events) {
            for (String treatment : treatments) {
                for (EvidenceDirection direction : EvidenceDirection.values()) {
                    result.addAll(reportHighestPerEventTreatmentDirection(evidences.stream()
                            .filter(x -> x.treatment().equals(treatment))
                            .filter(x -> x.direction().equals(direction))
                            .filter(x -> x.genomicEvent().equals(event))
                            .collect(Collectors.toList())));
                }
            }
        }

        return result;
    }

    @NotNull
    private static List<ProtectEvidence> reportHighestPerEventTreatmentDirection(@NotNull List<ProtectEvidence> evidences) {
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
            } else if (evidence.level() == highestOffLabel) {
                return highestOnLabel == null || evidence.level().isHigher(highestOnLabel);
            }
        }
        return false;
    }

    @Nullable
    @VisibleForTesting
    static EvidenceLevel highestReportableLevel(boolean isOnLabel, @NotNull List<ProtectEvidence> evidences) {
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
    public static List<ProtectEvidence> reportOnLabelTrialsOnly(@NotNull List<ProtectEvidence> evidences) {
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
        for (Knowledgebase source : evidence.sources()) {
            if (!TRIAL_SOURCES.contains(source)) {
                return false;
            }
        }
        return true;
    }
}
