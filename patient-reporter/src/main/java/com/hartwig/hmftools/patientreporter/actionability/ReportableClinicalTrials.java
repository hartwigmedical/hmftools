package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;

public final class ReportableClinicalTrials {

    private ReportableClinicalTrials() {
    }

    @NotNull
    public static List<EvidenceItem> reportableTrials(@NotNull List<EvidenceItem> evidenceItems) {
        Set<EvidenceItem> uniqueTrials = Sets.newHashSet();
        for (EvidenceItem evidence : evidenceItems) {
            if (evidence.source().isTrialSource()) {
                uniqueTrials.add(evidence);
            }
        }

        return Lists.newArrayList(uniqueTrials);
    }
}
