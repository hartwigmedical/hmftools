package com.hartwig.hmftools.orange.algo.protect;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

final class EvidenceSelector {

    static final Set<Knowledgebase> TRIAL_SOURCES = Sets.newHashSet(Knowledgebase.ACTIN, Knowledgebase.ICLUSION);

    private EvidenceSelector() {
    }

    @NotNull
    public static List<ProtectEvidence> reportableEvidence(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported() && !hasTrialSource(evidence.sources())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> reportableTrials(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported() && hasTrialSource(evidence.sources())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> unreportedEvidence(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (!evidence.reported() && !hasTrialSource(evidence.sources())) {
                filtered.add(evidence);
            }
        } return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> unreportedTrials(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (!evidence.reported() && hasTrialSource(evidence.sources())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    private static boolean hasTrialSource(@NotNull Iterable<ProtectSource> sources) {
        for (ProtectSource source : sources) {
            if (TRIAL_SOURCES.contains(source.name())) {
                return true;
            }
        }
        return false;
    }
}