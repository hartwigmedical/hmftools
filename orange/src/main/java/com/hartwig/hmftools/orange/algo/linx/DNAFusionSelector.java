package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.protect.EventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.algo.protect.EvidenceEvaluator;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class DNAFusionSelector {

    private DNAFusionSelector() {
    }

    @NotNull
    public static List<LinxFusion> selectInterestingUnreportedFusions(@NotNull List<LinxFusion> allFusions,
            @NotNull List<ProtectEvidence> evidences, @NotNull List<DriverGene> driverGenes) {
        List<LinxFusion> filtered = Lists.newArrayList();
        for (LinxFusion fusion : allFusions) {
            if (!fusion.reported()) {
                boolean hasEvidence = EvidenceEvaluator.hasEvidence(evidences, null, EventGenerator.fusionEvent(fusion));
                boolean hasReportedType = !fusion.reportedType().equals(KnownFusionType.NONE.toString());
                boolean isFusionOfOncogene = isInframeFusionWithOncogene(fusion, driverGenes);
                if (hasReportedType || hasEvidence || isFusionOfOncogene) {
                    filtered.add(fusion);
                }
            }
        }
        return filtered;
    }

    private static boolean isInframeFusionWithOncogene(@NotNull LinxFusion fusion, @NotNull List<DriverGene> driverGenes) {
        if (fusion.phased() != FusionPhasedType.INFRAME) {
            return false;
        }

        return isOncoDriverGene(driverGenes, fusion.geneStart()) || isOncoDriverGene(driverGenes, fusion.geneEnd());
    }

    private static boolean isOncoDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        DriverGene driver = findDriverGene(driverGenes, gene);
        return driver != null && driver.likelihoodType() == DriverCategory.ONCO;
    }

    @Nullable
    private static DriverGene findDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String geneToFind) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(geneToFind)) {
                return driverGene;
            }
        }
        return null;
    }
}
