package com.hartwig.hmftools.patientreporter.actionability;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.EvidenceItem;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportableEvidenceItems {
    private static final Logger LOGGER = LogManager.getLogger(ReportableEvidenceItems.class);

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

        HashSet<EvidenceItem> filteringReportableItems = new HashSet<>(reportableItems);
        ArrayList<EvidenceItem> filteredReportableItems = new ArrayList<>(filteringReportableItems);

        List<EvidenceItem> evidenceUnique = Lists.newArrayList();
        for (EvidenceItem distinctTrials : filteredReportableItems) {
            EvidenceItem valueEvidenceUnique = evidenceUnique.iterator().next();
            if (evidenceUnique.isEmpty()) {
                evidenceUnique.add(distinctTrials);
            } else if (distinctTrials.event().equals(valueEvidenceUnique.event()) && distinctTrials.drug()
                    .equals(valueEvidenceUnique.drug()) && distinctTrials.response()
                    .equals(valueEvidenceUnique.response()) && !valueEvidenceUnique.isOnLabel()) {
                evidenceUnique.remove(distinctTrials);
            } else if (distinctTrials.event().equals(valueEvidenceUnique.event()) && distinctTrials.drug()
                    .equals(valueEvidenceUnique.drug()) && distinctTrials.response()
                    .equals(valueEvidenceUnique.response()) && !distinctTrials.isOnLabel()) {
                evidenceUnique.remove(distinctTrials);
            } else {
                evidenceUnique.add(distinctTrials);
            }
        }
        return evidenceUnique;
    }
}