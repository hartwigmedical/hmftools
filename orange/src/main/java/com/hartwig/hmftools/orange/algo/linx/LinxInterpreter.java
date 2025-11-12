package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

public class LinxInterpreter
{
    private final EnsemblDataCache mEnsemblDataCache;

    public LinxInterpreter(final EnsemblDataCache ensemblDataCache)
    {
        mEnsemblDataCache = ensemblDataCache;
    }

    public LinxRecord interpret(final LinxData linx)
    {
        LOGGER.info("Analysing linx data");

        LinxBreakendInterpreter somaticBreakendInterpreter = new LinxBreakendInterpreter(
                linx.allSomaticSvAnnotations(),
                mEnsemblDataCache);

        LinxBreakendInterpreter germlineBreakendInterpreter = new LinxBreakendInterpreter(
                Objects.requireNonNullElse(linx.allGermlineSvAnnotations(), List.of()),
                mEnsemblDataCache);

        return ImmutableLinxRecord.builder()
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDrivers(), LinxConversion::convert))
                .allSomaticStructuralVariants(ConversionUtil.mapToIterable(linx.allSomaticSvAnnotations(), LinxConversion::convert))
                .allSomaticFusions(ConversionUtil.mapToIterable(linx.allSomaticFusions(), LinxConversion::convert))
                .reportableSomaticFusions(ConversionUtil.mapToIterable(linx.reportableSomaticFusions(), LinxConversion::convert))
                .allSomaticBreakends(ConversionUtil.mapToIterable(linx.allSomaticBreakends(), somaticBreakendInterpreter::interpret))
                .reportableSomaticBreakends(ConversionUtil.mapToIterable(linx.reportableSomaticBreakends(), somaticBreakendInterpreter::interpret))
                .somaticHomozygousDisruptions(ConversionUtil.mapToIterable(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .allGermlineStructuralVariants(ConversionUtil.mapToIterable(linx.allGermlineSvAnnotations(), LinxConversion::convert))
                .allGermlineBreakends(ConversionUtil.mapToIterable(linx.allGermlineBreakends(), germlineBreakendInterpreter::interpret))
                .reportableGermlineBreakends(ConversionUtil.mapToIterable(linx.reportableGermlineBreakends(), germlineBreakendInterpreter::interpret))
                .germlineHomozygousDisruptions(ConversionUtil.mapToIterable(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }
}
