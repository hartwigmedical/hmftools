package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectSource;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrialFactory {

    private ClinicalTrialFactory() {
    }

    @NotNull
    public static List<ProtectEvidence> extractOnLabelTrials(@NotNull List<ProtectEvidence> evidenceItems) {
        List<ProtectEvidence> trials = Lists.newArrayList();
        for (ProtectEvidence evidence : evidenceItems) {
            for (ProtectSource protectSource: evidence.protectSources()) {
                if (protectSource.source() == Knowledgebase.ICLUSION && evidence.onLabel()) {
                    trials.add(evidence);
                }
            }

        }
        return trials;
    }
}