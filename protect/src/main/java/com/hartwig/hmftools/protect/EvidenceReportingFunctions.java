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
    private static final EvidenceLevel MAX_REPORTABLE_LEVEL = EvidenceLevel.B;

    private EvidenceReportingFunctions() {
    }

    @NotNull
    public static List<ProtectEvidence> reportHighestLevelEvidence(@NotNull List<ProtectEvidence> evidence) {
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
        if (evidence.reported() && evidence.level().ordinal() <= MAX_REPORTABLE_LEVEL.ordinal()) {
            if (evidence.onLabel()) {
                assert highestOnLabel != null;
                return evidence.level() == highestOnLabel;
            } else if (evidence.level() == highestOffLabel){
                return highestOnLabel == null || evidence.level().ordinal() < highestOnLabel.ordinal();
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
                if (highest == null || evidence.level().ordinal() < highest.ordinal()) {
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
            if (!hasAtLeastOneNoneTrialSource(evidence)) {
                result.add(ImmutableProtectEvidence.builder().from(evidence).reported(evidence.reported() && evidence.onLabel()).build());
            } else {
                result.add(evidence);
            }
        }
        return result;
    }

    private static boolean hasAtLeastOneNoneTrialSource(@NotNull ProtectEvidence evidence) {
        for (Knowledgebase source : evidence.sources()) {
            if (!TRIAL_SOURCES.contains(source)) {
                return true;
            }
        }
        return false;
    }
}
