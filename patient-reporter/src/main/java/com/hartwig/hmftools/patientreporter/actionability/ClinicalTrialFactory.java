package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrialFactory {

    private ClinicalTrialFactory() {
    }

    @NotNull
    public static List<ClinicalTrial> extractOnLabelTrials(@NotNull List<EvidenceItem> evidenceItems) {
        List<ClinicalTrial> trials = Lists.newArrayList();
        for (EvidenceItem evidence : evidenceItems) {
            if (evidence.source().isTrialSource() && evidence.isOnLabel()) {
                trials.add(toClinicalTrial(evidence));
            }
        }

        return trials;
    }

    @NotNull
    private static ClinicalTrial toClinicalTrial(@NotNull EvidenceItem evidenceItem) {
        return ImmutableClinicalTrial.builder()
                .event(evidenceItem.event())
                .acronym(evidenceItem.drug())
                .source(evidenceItem.source())
                .reference(evidenceItem.reference())
                .isOnLabel(evidenceItem.isOnLabel())
                .cancerType(evidenceItem.cancerType())
                .build();
    }
}
