package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CheckEvidenceCnv {

    private static final Logger LOGGER = LogManager.getLogger(CheckEvidenceCnv.class);

    private CheckEvidenceCnv() {
    }

    static Map<GeneCopyNumber, List<EvidenceItem>> checkingAndFilterForEvidenceInDriverCatalog(@NotNull List<ReportableGainLoss> reportableGainLosses,
            Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber, @NotNull List<EvidenceItem> allEvidenceForCopyNumbers) {
        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        Map<GeneCopyNumber, List<EvidenceItem>> filterEvidenceMap = Maps.newHashMap();

        // remove evidence for not reportable CNV
        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), allEvidenceForCopyNumbers) && !reportableGenes.contains(geneCopyNumber.gene())) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            } else {
                filterEvidenceMap.put(entry.getKey(), evidencePerGeneCopyNumber.get(geneCopyNumber));
            }
        }

        return filterEvidenceMap;
    }
}