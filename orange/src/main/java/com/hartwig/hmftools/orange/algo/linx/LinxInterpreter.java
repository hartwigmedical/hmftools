package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

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
                .allSomaticStructuralVariants(ConversionUtil.convertCollection(linx.allSomaticStructuralVariants(), LinxConversion::convert))
                .allSomaticFusions(ConversionUtil.convertCollection(linx.allSomaticFusions(), LinxConversion::convert))
                .reportableSomaticFusions(ConversionUtil.convertCollection(linx.reportableSomaticFusions(), LinxConversion::convert))
                .additionalSuspectSomaticFusions(ConversionUtil.convertCollection(additionalSuspectSomaticFusions, LinxConversion::convert))
                .allSomaticBreakends(ConversionUtil.convertCollection(linx.allSomaticBreakends(), LinxConversion::convert))
                .reportableSomaticBreakends(ConversionUtil.convertCollection(linx.reportableSomaticBreakends(), LinxConversion::convert))
                .additionalSuspectSomaticBreakends(ConversionUtil.convertCollection(additionalSuspectSomaticBreakends, LinxConversion::convert))
                .somaticHomozygousDisruptions(ConversionUtil.convertCollection(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .allGermlineStructuralVariants(ConversionUtil.convertCollection(linx.allGermlineStructuralVariants(), LinxConversion::convert))
                .allGermlineBreakends(ConversionUtil.convertCollection(linx.allGermlineBreakends(), LinxConversion::convert))
                .reportableGermlineBreakends(ConversionUtil.convertCollection(linx.reportableGermlineBreakends(), LinxConversion::convert))
                .germlineHomozygousDisruptions(ConversionUtil.convertCollection(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }
}
