package com.hartwig.hmftools.common.actionability.cnv;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.reportablegenomicalterations.ReportableGainLoss;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CheckEvidenceCnv {
    private static final Logger LOGGER = LogManager.getLogger(CheckEvidenceCnv.class);

    private CheckEvidenceCnv() {

    }

    public static void checkingForEvidenceInDriverCatalog(@NotNull List<ReportableGainLoss> reportableGainLosses,
            Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber, @NotNull List<EvidenceItem> allEvidenceForCopyNumbers) {
        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), allEvidenceForCopyNumbers) && !reportableGenes.contains(geneCopyNumber.gene())) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            }
        }
    }
}