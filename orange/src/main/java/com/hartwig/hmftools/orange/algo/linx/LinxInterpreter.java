package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.orange.conversion.LinxConversion;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Objects;

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
    public LinxRecord interpret(@NotNull LinxData linx) {
        LOGGER.info("Analysing linx data");
        List<com.hartwig.hmftools.common.linx.LinxFusion> additionalSuspectSomaticFusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(linx.allSomaticFusions(), driverGenes);
        LOGGER.info(" Found an additional {} suspect somatic fusions that are potentially interesting",
                additionalSuspectSomaticFusions.size());

        List<com.hartwig.hmftools.common.linx.LinxBreakend> additionalSuspectSomaticBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allSomaticBreakends(),
                        linx.reportableSomaticFusions(),
                        knownFusionCache);
        LOGGER.info(" Found an additional {} suspect somatic breakends that are potentially interesting",
                additionalSuspectSomaticBreakends.size());

        return ImmutableLinxRecord.builder()
                .allSomaticStructuralVariants(() -> linx.allSomaticStructuralVariants().stream().map(LinxConversion::convert).iterator())
                .allGermlineStructuralVariants(() -> argumentOrEmpty(linx.allGermlineStructuralVariants()).stream().map(LinxConversion::convert).iterator())
                .allSomaticFusions(() -> linx.allSomaticFusions().stream().map(LinxConversion::convert).iterator())
                .reportableSomaticFusions(() -> linx.reportableSomaticFusions().stream().map(LinxConversion::convert).iterator())
                .additionalSuspectSomaticFusions(() -> additionalSuspectSomaticFusions.stream().map(LinxConversion::convert).iterator())
                .allSomaticBreakends(() -> linx.allSomaticBreakends().stream().map(LinxConversion::convert).iterator())
                .allGermlineBreakends(() -> argumentOrEmpty(linx.allGermlineBreakends()).stream().map(LinxConversion::convert).iterator())
                .reportableSomaticBreakends(() -> linx.reportableSomaticBreakends().stream().map(LinxConversion::convert).iterator())
                .reportableGermlineBreakends(() -> argumentOrEmpty(linx.reportableGermlineBreakends()).stream().map(LinxConversion::convert).iterator())
                .additionalSuspectSomaticBreakends(() -> additionalSuspectSomaticBreakends.stream().map(LinxConversion::convert).iterator())
                .somaticHomozygousDisruptions(() -> linx.somaticHomozygousDisruptions().stream().map(LinxConversion::convert).iterator())
                .build();
    }

    private static <T> List<T> argumentOrEmpty(List<T> obj) {
        return Objects.requireNonNullElseGet(obj, List::of);
    }
}
