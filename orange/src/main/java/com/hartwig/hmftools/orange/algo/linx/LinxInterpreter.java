package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LinxInterpreter
{
    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final List<StructuralVariant> allSomaticStructuralVariants;
    @Nullable
    private final List<StructuralVariant> allGermlineStructuralVariants;
    @NotNull
    private final List<StructuralVariant> allInferredSomaticStructuralVariants;
    @Nullable
    private final List<StructuralVariant> allInferredGermlineStructuralVariants;
    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public LinxInterpreter(
            @NotNull final List<DriverGene> driverGenes,
            @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final List<StructuralVariant> allSomaticStructuralVariants,
            @Nullable final List<StructuralVariant> allGermlineStructuralVariants,
            @NotNull final List<StructuralVariant> allInferredSomaticStructuralVariants,
            @Nullable final List<StructuralVariant> allInferredGermlineStructuralVariants,
            @NotNull final EnsemblDataCache ensemblDataCache
    )
    {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.allSomaticStructuralVariants = allSomaticStructuralVariants;
        this.allGermlineStructuralVariants = allGermlineStructuralVariants;
        this.allInferredSomaticStructuralVariants = allInferredSomaticStructuralVariants;
        this.allInferredGermlineStructuralVariants = allInferredGermlineStructuralVariants;
        this.ensemblDataCache = ensemblDataCache;
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

        LinxBreakendInterpreter somaticBreakendInterpreter = new LinxBreakendInterpreter(
                combine(allSomaticStructuralVariants, allInferredSomaticStructuralVariants),
                linx.allSomaticStructuralVariants(),
                ensemblDataCache);

        LinxBreakendInterpreter germlineBreakendInterpreter = new LinxBreakendInterpreter(
                combine(allGermlineStructuralVariants, allInferredGermlineStructuralVariants),
                Objects.requireNonNullElse(linx.allGermlineStructuralVariants(), List.of()),
                ensemblDataCache);

        return ImmutableLinxRecord.builder()
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDrivers(), LinxConversion::convert))
                .allSomaticStructuralVariants(ConversionUtil.mapToIterable(linx.allSomaticStructuralVariants(), LinxConversion::convert))
                .allSomaticFusions(ConversionUtil.mapToIterable(linx.allSomaticFusions(), LinxConversion::convert))
                .reportableSomaticFusions(ConversionUtil.mapToIterable(linx.reportableSomaticFusions(), LinxConversion::convert))
                .additionalSuspectSomaticFusions(ConversionUtil.mapToIterable(additionalSuspectSomaticFusions, LinxConversion::convert))
                .additionalViableSomaticFusions(ConversionUtil.mapToIterable(additionalViableSomaticFusions, LinxConversion::convert))
                .allSomaticBreakends(ConversionUtil.mapToIterable(linx.allSomaticBreakends(), somaticBreakendInterpreter::interpret))
                .reportableSomaticBreakends(ConversionUtil.mapToIterable(linx.reportableSomaticBreakends(), somaticBreakendInterpreter::interpret))
                .additionalSuspectSomaticBreakends(ConversionUtil.mapToIterable(additionalSuspectSomaticBreakends, somaticBreakendInterpreter::interpret))
                .somaticHomozygousDisruptions(ConversionUtil.mapToIterable(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .allGermlineStructuralVariants(ConversionUtil.mapToIterable(linx.allGermlineStructuralVariants(), LinxConversion::convert))
                .allGermlineBreakends(ConversionUtil.mapToIterable(linx.allGermlineBreakends(), germlineBreakendInterpreter::interpret))
                .reportableGermlineBreakends(ConversionUtil.mapToIterable(linx.reportableGermlineBreakends(), germlineBreakendInterpreter::interpret))
                .germlineHomozygousDisruptions(ConversionUtil.mapToIterable(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }

    @NotNull
    private List<StructuralVariant> combine(@Nullable List<StructuralVariant> svList1, @Nullable List<StructuralVariant> svList2)
    {
        return Stream.concat(
                        Optional.ofNullable(svList1).orElse(List.of()).stream(),
                        Optional.ofNullable(svList2).orElse(List.of()).stream())
                .collect(Collectors.toList());
    }
}
