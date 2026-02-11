package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

public class LinxInterpreter
{
    private final CytoBands mCytoBands;

    public LinxInterpreter(final CytoBands cytoBands)
    {
        mCytoBands = cytoBands;
    }

    public LinxRecord interpret(final LinxData linx)
    {
        LOGGER.info("Analysing Linx data");

        LinxBreakendInterpreter somaticBreakendInterpreter = new LinxBreakendInterpreter(linx.somaticSvAnnotations(), mCytoBands);

        LinxBreakendInterpreter germlineBreakendInterpreter = new LinxBreakendInterpreter(
                Objects.requireNonNullElse(linx.germlineSvAnnotations(), List.of()), mCytoBands);

        return ImmutableLinxRecord.builder()
                .somaticStructuralVariants(ConversionUtil.mapToIterable(linx.somaticSvAnnotations(), LinxConversion::convert))
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDrivers(), LinxConversion::convert))
                .fusions(ConversionUtil.mapToIterable(linx.fusions(), LinxConversion::convert))
                .somaticBreakends(ConversionUtil.mapToIterable(linx.somaticBreakends(), somaticBreakendInterpreter::interpret))
                .somaticHomozygousDisruptions(ConversionUtil.mapToIterable(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .germlineStructuralVariants(ConversionUtil.mapToIterable(linx.germlineSvAnnotations(), LinxConversion::convert))
                .germlineBreakends(ConversionUtil.mapToIterable(linx.germlineBreakends(), germlineBreakendInterpreter::interpret))
                .germlineHomozygousDisruptions(ConversionUtil.mapToIterable(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }
}
