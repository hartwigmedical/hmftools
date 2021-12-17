package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class EvidenceSelector {

    private EvidenceSelector() {
    }

    public static boolean hasEvidence(@NotNull List<ProtectEvidence> evidences, @NotNull String genomicEvent) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.genomicEvent().equals(genomicEvent)) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    public static List<ProtectEvidence> reported(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.reported()) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> unreported(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (!evidence.reported()) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> iclusionOnly(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (evidence.sources().contains(Knowledgebase.ICLUSION)) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static List<ProtectEvidence> noIclusion(@NotNull List<ProtectEvidence> evidences) {
        List<ProtectEvidence> filtered = Lists.newArrayList();
        for (ProtectEvidence evidence : evidences) {
            if (!evidence.sources().contains(Knowledgebase.ICLUSION)) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    public static Map<String, List<ProtectEvidence>> buildTreatmentMap(@NotNull List<ProtectEvidence> evidences, boolean requireOnLabel) {
        Map<String, List<ProtectEvidence>> evidencePerTreatmentMap = Maps.newHashMap();

        for (ProtectEvidence evidence : evidences) {
            if (evidence.onLabel() == requireOnLabel) {
                List<ProtectEvidence> treatmentEvidences = evidencePerTreatmentMap.get(evidence.treatment());
                if (treatmentEvidences == null) {
                    treatmentEvidences = Lists.newArrayList();
                }
                if (!hasHigherOrEqualEvidenceForEventAndTreatment(treatmentEvidences, evidence)) {
                    treatmentEvidences.add(evidence);
                }
                evidencePerTreatmentMap.put(evidence.treatment(), treatmentEvidences);
            }
        }
        return evidencePerTreatmentMap;
    }

    private static boolean hasHigherOrEqualEvidenceForEventAndTreatment(@NotNull List<ProtectEvidence> evidences,
            @NotNull ProtectEvidence evidenceToCheck) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(evidenceToCheck.treatment()) && evidence.genomicEvent()
                    .equals(evidenceToCheck.genomicEvent())) {
                if (!evidenceToCheck.level().isHigher(evidence.level())) {
                    return true;
                }
            }
        }
        return false;
    }
}
