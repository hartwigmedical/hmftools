package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.orange.algo.selection.FusionSelector;

import org.jetbrains.annotations.NotNull;

public final class LinxInterpreter {

    private LinxInterpreter() {
    }

    @NotNull
    public static LinxInterpretedData interpret(@NotNull LinxData linx, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<DriverGene> driverGenes) {
        return ImmutableLinxInterpretedData.builder()
                .allFusions(linx.allFusions())
                .reportableFusions(linx.reportableFusions())
                .potentiallyInterestingFusions(FusionSelector.selectInterestingUnreportedFusions(linx.allFusions(),
                        evidences,
                        driverGenes))
                .geneDisruptions(linx.geneDisruptions())
                .homozygousDisruptions(linx.homozygousDisruptions())
                .drivers(linx.drivers())
                .allGermlineDisruptions(linx.allGermlineDisruptions())
                .reportableGermlineDisruptions(linx.reportableGermlineDisruptions())
                .build();
    }
}
