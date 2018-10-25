package com.hartwig.hmftools.patientreporter.actionability;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.jetbrains.annotations.NotNull;

public final class ReportableEvidenceItems {

    private ReportableEvidenceItems() {
    }

    @NotNull
    public static List<EvidenceItem> reportableEvidenceItems(@NotNull List<EvidenceItem> evidenceItems) {
        List<EvidenceItem> reportableItems = Lists.newArrayList();
        for (EvidenceItem evidence : evidenceItems) {
            if (!evidence.source().isTrialSource()) {
                reportableItems.add(evidence);
            }
        }
        return reportableItems;
    }
}
