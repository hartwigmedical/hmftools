package com.hartwig.hmftools.common.actionability;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class ReportableEvidenceItem {

    @NotNull
    public static List<EvidenceItem> extractAllReportableEvidenceItems(@NotNull List<EvidenceItem> allEvidenceItems) {
        List<EvidenceItem> allEvidenceItemsFiltered = Lists.newArrayList();
        for (EvidenceItem evidenceItem: allEvidenceItems) {
            if (!evidenceItem.event().split(" ")[0].equals("TP53") && !evidenceItem.drug().equals("Tamoxifen")) {
                allEvidenceItemsFiltered.add(ImmutableEvidenceItem.builder()
                        .event(evidenceItem.event())
                        .isOnLabel(evidenceItem.isOnLabel())
                        .reference(evidenceItem.reference())
                        .source(evidenceItem.source())
                        .drug(evidenceItem.drug())
                        .drugsType(evidenceItem.drugsType())
                        .level(evidenceItem.level())
                        .response(evidenceItem.response())
                        .cancerType(evidenceItem.cancerType())
                        .scope(evidenceItem.scope()).build());
            }
        }
        return allEvidenceItemsFiltered;
    }
}
