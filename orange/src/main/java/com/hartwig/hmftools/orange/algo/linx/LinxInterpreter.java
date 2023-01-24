package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxFusion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LinxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LinxInterpreter.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public LinxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache) {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
    }

    @NotNull
    public LinxInterpretedData interpret(@NotNull LinxData linx) {
        List<LinxFusion> additionalSuspectSomaticFusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(linx.allSomaticFusions(), driverGenes);
        LOGGER.info(" Found an additional {} suspect fusions that are potentially interesting", additionalSuspectSomaticFusions.size());

        List<LinxBreakend> additionalSuspectSomaticBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allSomaticBreakends(),
                        linx.reportableSomaticFusions(),
                        knownFusionCache);
        LOGGER.info(" Found an additional {} suspect breakends that are potentially interesting", additionalSuspectSomaticBreakends.size());

        return ImmutableLinxInterpretedData.builder()
                .allSomaticStructuralVariants(linx.allSomaticStructuralVariants())
                .allSomaticFusions(linx.allSomaticFusions())
                .reportableSomaticFusions(linx.reportableSomaticFusions())
                .additionalSuspectSomaticFusions(additionalSuspectSomaticFusions)
                .allSomaticBreakends(linx.allSomaticBreakends())
                .reportableSomaticBreakends(linx.reportableSomaticBreakends())
                .additionalSuspectSomaticBreakends(additionalSuspectSomaticBreakends)
                .somaticHomozygousDisruptions(linx.somaticHomozygousDisruptions())
                .allGermlineStructuralVariants(linx.allGermlineStructuralVariants())
                .allGermlineBreakends(linx.allGermlineBreakends())
                .reportableGermlineBreakends(linx.reportableGermlineBreakends())
                .allGermlineDisruptions(linx.allGermlineDisruptions())
                .reportableGermlineDisruptions(linx.reportableGermlineDisruptions())
                .build();
    }
}
