package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CheckEvidenceCnv {

    private static final Logger LOGGER = LogManager.getLogger(CheckEvidenceCnv.class);

    private CheckEvidenceCnv() {
    }

    static Map<GeneCopyNumber, List<EvidenceItem>> checkingAndFilterForEvidenceInDriverCatalog(
            @NotNull List<ReportableGainLoss> reportableGainLosses, Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber) {
        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        Map<GeneCopyNumber, List<EvidenceItem>> filterEvidenceMap = Maps.newHashMap();
        Map<GeneCopyNumber, List<EvidenceItem>> evidenceMapNonReportable = Maps.newHashMap();

        // remove evidence for not reportable CNV
        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (reportableGenes.contains(geneCopyNumber.gene())) {
                filterEvidenceMap.put(entry.getKey(), evidencePerGeneCopyNumber.get(geneCopyNumber));
            } else {
                evidenceMapNonReportable.put(entry.getKey(), evidencePerGeneCopyNumber.get(geneCopyNumber));

            }
        }

        List<EvidenceItem> evidenceNonReportableItem = ReportableEvidenceItemFactory.toList(evidenceMapNonReportable);
        Set<String> uniqueEventsNonReportableEvidence = Sets.newHashSet();
        for (EvidenceItem item : evidenceNonReportableItem) {
            if (!item.level().equals(EvidenceLevel.LEVEL_A) || !item.level().equals(EvidenceLevel.LEVEL_B)) {
                uniqueEventsNonReportableEvidence.add(item.event());
            }
        }

        for (String event: uniqueEventsNonReportableEvidence ) {
            LOGGER.warn("Copy evidence not reported for event " + event + "!");
        }
        return filterEvidenceMap;
    }
}