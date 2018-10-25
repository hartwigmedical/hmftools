package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;

public final class ReportableClinicalTrials {

    private ReportableClinicalTrials() {
    }

    @NotNull
    public static List<EvidenceItem> reportableTrials(@NotNull List<EvidenceItem> evidenceItems) {
        List<EvidenceItem> trials = Lists.newArrayList();
        for (EvidenceItem evidence : evidenceItems) {
            if (evidence.source().isTrialSource()) {
                trials.add(evidence);
            }
        }
        return trials;
    }
}
