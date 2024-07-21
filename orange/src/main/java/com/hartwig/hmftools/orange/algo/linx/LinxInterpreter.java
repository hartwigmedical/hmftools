package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.jetbrains.annotations.NotNull;

public class LinxInterpreter
{
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public LinxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache)
    {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
    }

    @NotNull
    public LinxRecord interpret(@NotNull LinxData linx)
    {
        LOGGER.info("Analysing linx data");
        List<LinxFusion> additionalSuspectSomaticFusions =
                DnaFusionSelector.selectInterestingUnreportedFusions(linx.allSomaticFusions(), driverGenes);
        LOGGER.info(" Found an additional {} suspect somatic fusions that are potentially interesting",
                additionalSuspectSomaticFusions.size());

        List<LinxFusion> additionalViableSomaticFusions =
                DnaFusionSelector.selectAdditionalViableSomaticFusions(linx.allSomaticFusions(), additionalSuspectSomaticFusions);
        LOGGER.info(" Found an additional {} viable somatic fusions", additionalViableSomaticFusions.size());

        List<LinxBreakend> additionalSuspectSomaticBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allSomaticBreakends(),
                        linx.reportableSomaticFusions(),
                        knownFusionCache);
        LOGGER.info(" Found an additional {} suspect somatic breakends that are potentially interesting",
                additionalSuspectSomaticBreakends.size());

        return ImmutableLinxRecord.builder()
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDrivers(), LinxConversion::convert))
                .allSomaticStructuralVariants(ConversionUtil.mapToIterable(linx.allSomaticStructuralVariants(), LinxConversion::convert))
                .allSomaticFusions(ConversionUtil.mapToIterable(linx.allSomaticFusions(), LinxConversion::convert))
                .reportableSomaticFusions(ConversionUtil.mapToIterable(linx.reportableSomaticFusions(), LinxConversion::convert))
                .additionalSuspectSomaticFusions(ConversionUtil.mapToIterable(additionalSuspectSomaticFusions, LinxConversion::convert))
                .additionalViableSomaticFusions(ConversionUtil.mapToIterable(additionalViableSomaticFusions, LinxConversion::convert))
                .allSomaticBreakends(ConversionUtil.mapToIterable(linx.allSomaticBreakends(), LinxConversion::convert))
                .reportableSomaticBreakends(ConversionUtil.mapToIterable(linx.reportableSomaticBreakends(), LinxConversion::convert))
                .additionalSuspectSomaticBreakends(ConversionUtil.mapToIterable(additionalSuspectSomaticBreakends, LinxConversion::convert))
                .somaticHomozygousDisruptions(ConversionUtil.mapToIterable(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .allGermlineStructuralVariants(ConversionUtil.mapToIterable(linx.allGermlineStructuralVariants(), LinxConversion::convert))
                .allGermlineBreakends(ConversionUtil.mapToIterable(linx.allGermlineBreakends(), LinxConversion::convert))
                .reportableGermlineBreakends(ConversionUtil.mapToIterable(linx.reportableGermlineBreakends(), LinxConversion::convert))
                .germlineHomozygousDisruptions(ConversionUtil.mapToIterable(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }
}
