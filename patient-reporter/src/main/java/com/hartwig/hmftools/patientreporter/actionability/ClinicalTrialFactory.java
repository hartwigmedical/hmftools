package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableClinicalTrial;

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
                .scope(evidenceItem.scope())
                .type(evidenceItem.type())
                .gene(evidenceItem.gene())
                .chromosome(evidenceItem.chromosome())
                .position(evidenceItem.position())
                .ref(evidenceItem.ref())
                .alt(evidenceItem.alt())
                .cnvType(evidenceItem.cnvType())
                .fusionFiveGene(evidenceItem.fusionFiveGene())
                .fusionThreeGene(evidenceItem.fusionThreeGene())
                .build();
    }
}
