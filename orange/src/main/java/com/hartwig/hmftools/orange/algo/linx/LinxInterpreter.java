package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.GeneDisruption;
import com.hartwig.hmftools.common.linx.GeneDisruptionFactory;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LinxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LinxInterpreter.class);

    private LinxInterpreter() {
    }

    @NotNull
    public static LinxInterpretedData interpret(@NotNull LinxData linx, @NotNull List<ProtectEvidence> evidences,
            @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache) {
        List<LinxFusion> additionalSuspectFusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(linx.allFusions(), evidences, driverGenes);
        LOGGER.info(" Found an additional {} suspect fusions that are potentially interesting", additionalSuspectFusions.size());

        List<LinxBreakend> additionalSuspectBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allBreakends(), linx.reportableFusions(), knownFusionCache);
        List<GeneDisruption> additionalSuspectDisruptions =
                GeneDisruptionFactory.convert(additionalSuspectBreakends, linx.allStructuralVariants());
        LOGGER.info(" Found an additional {} suspect disruptions that are potentially interesting", additionalSuspectDisruptions.size());

        return ImmutableLinxInterpretedData.builder()
                .allStructuralVariants(linx.allStructuralVariants())
                .allFusions(linx.allFusions())
                .reportableFusions(linx.reportableFusions())
                .additionalSuspectFusions(additionalSuspectFusions)
                .allBreakends(linx.allBreakends())
                .reportableGeneDisruptions(linx.reportableGeneDisruptions())
                .additionalSuspectDisruptions(additionalSuspectDisruptions)
                .homozygousDisruptions(linx.homozygousDisruptions())
                .allGermlineDisruptions(linx.allGermlineDisruptions())
                .reportableGermlineDisruptions(linx.reportableGermlineDisruptions())
                .build();
    }
}
