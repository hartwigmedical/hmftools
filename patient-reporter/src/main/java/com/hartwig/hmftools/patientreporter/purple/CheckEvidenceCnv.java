package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CheckEvidenceCnv {

    private static final Logger LOGGER = LogManager.getLogger(CheckEvidenceCnv.class);

    private CheckEvidenceCnv() {
    }

    @NotNull
    static Map<ReportableGainLoss, List<EvidenceItem>> checkAndFilterForEvidenceInDriverCatalog(
            @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull Map<ReportableGainLoss, List<EvidenceItem>> evidencePerGeneCopyNumber) {
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        Map<ReportableGainLoss, List<EvidenceItem>> filteredEvidenceMap = Maps.newHashMap();
        Map<ReportableGainLoss, List<EvidenceItem>> evidenceMapNonReportable = Maps.newHashMap();

        // Remove evidence for not reportable CNV
        for (Map.Entry<ReportableGainLoss, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            ReportableGainLoss geneCopyNumber = entry.getKey();
            if (reportableGenes.contains(geneCopyNumber.gene())) {
                filteredEvidenceMap.put(entry.getKey(), evidencePerGeneCopyNumber.get(geneCopyNumber));
            } else {
                evidenceMapNonReportable.put(entry.getKey(), evidencePerGeneCopyNumber.get(geneCopyNumber));
            }
        }

        // Report a warning for all events that are filtered out with A or B level evidence.
        Set<String> uniqueEventsNonReportableEvidence = Sets.newHashSet();
        for (Map.Entry<ReportableGainLoss, List<EvidenceItem>> entry : evidenceMapNonReportable.entrySet()) {
            for (EvidenceItem item : entry.getValue()) {
                if (item.level() == EvidenceLevel.LEVEL_A || item.level() == EvidenceLevel.LEVEL_B) {
                    uniqueEventsNonReportableEvidence.add(item.event());
                }
            }
        }

        for (String event : uniqueEventsNonReportableEvidence) {
            LOGGER.warn("Copy evidence not reported for event {}!", event);
        }

        return filteredEvidenceMap;
    }
}