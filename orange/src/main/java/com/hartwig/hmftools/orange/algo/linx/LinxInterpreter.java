package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.algo.selection.FusionSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LinxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LinxInterpreter.class);

    private LinxInterpreter() {
    }

    @NotNull
    public static LinxInterpretedData interpret(@NotNull LinxData linx, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<DriverGene> driverGenes) {
        List<LinxFusion> additionalSuspectFusions =
                FusionSelector.selectInterestingUnreportedFusions(linx.allFusions(), evidences, driverGenes);
        LOGGER.info(" Found an additional {} suspect fusions that are potentially interesting", additionalSuspectFusions.size());

        return ImmutableLinxInterpretedData.builder()
                .allFusions(linx.allFusions())
                .reportableFusions(linx.reportableFusions())
                .additionalSuspectFusions(additionalSuspectFusions)
                .geneDisruptions(linx.reportableGeneDisruptions())
                .homozygousDisruptions(linx.homozygousDisruptions())
                .drivers(linx.drivers())
                .allGermlineDisruptions(linx.allGermlineDisruptions())
                .reportableGermlineDisruptions(linx.reportableGermlineDisruptions())
                .build();
    }
}
