package com.hartwig.hmftools.patientreporter.actionability;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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
        return new ArrayList<>(filteringReportableItems);
    }
}
