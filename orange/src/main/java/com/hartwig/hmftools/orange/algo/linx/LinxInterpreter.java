package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.linx.DisruptionFactory.createDisruptions;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.orange.algo.util.FindingKeys;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

public class LinxInterpreter
{
    private final EnsemblDataCache mEnsemblDataCache;

    public LinxInterpreter(final EnsemblDataCache ensemblDataCache)
    {
        mEnsemblDataCache = ensemblDataCache;
    }

    public LinxRecord interpret(final LinxData linx, final PurpleQC purpleQC)
    {
        LOGGER.info("Analysing linx data");

        boolean hasReliablePurity = !purpleQC.status().contains(PurpleQCStatus.FAIL_NO_TUMOR);

        LinxBreakendInterpreter somaticBreakendInterpreter = new LinxBreakendInterpreter(
                linx.allSomaticSvAnnotations(),
                mEnsemblDataCache);

        List<LinxSvAnnotation> allSomaticSvAnnotations = ConversionUtil.mapToList(linx.allSomaticSvAnnotations(), LinxConversion::convert);
        List<LinxBreakend> driverSomaticBreakends = ConversionUtil.mapToList(linx.driverSomaticBreakends(), somaticBreakendInterpreter::interpret);
        List<Disruption> driverSomaticDisruptions = createDisruptions(FindingKeys.SampleType.SOMATIC, driverSomaticBreakends,
                allSomaticSvAnnotations, hasReliablePurity);

        List<LinxSvAnnotation> allGermlineSvAnnotations = null;
        List<LinxBreakend> driverGermlineBreakends = null;
        List<LinxBreakend> otherGermlineBreakends = null;
        List<Disruption> driverGermlineDisruptions = null;

        if(linx.allGermlineSvAnnotations() != null)
        {
            LinxBreakendInterpreter germlineBreakendInterpreter = new LinxBreakendInterpreter(
                    Objects.requireNonNull(linx.allGermlineSvAnnotations()), mEnsemblDataCache);

            allGermlineSvAnnotations = ConversionUtil.mapToList(linx.allGermlineSvAnnotations(), LinxConversion::convert);
            driverGermlineBreakends = ConversionUtil.mapToList(linx.driverGermlineBreakends(), germlineBreakendInterpreter::interpret);
            otherGermlineBreakends = ConversionUtil.mapToList(linx.allGermlineBreakends(), germlineBreakendInterpreter::interpret);
            driverGermlineDisruptions = createDisruptions(FindingKeys.SampleType.GERMLINE, driverGermlineBreakends,
                    allGermlineSvAnnotations, hasReliablePurity);
        }

        return ImmutableLinxRecord.builder()
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDrivers(), LinxConversion::convert))
                .allSomaticStructuralVariants(allSomaticSvAnnotations)
                .allSomaticFusions(ConversionUtil.mapToIterable(linx.allSomaticFusions(), LinxConversion::convert))
                .otherSomaticBreakends(ConversionUtil.mapToIterable(linx.allSomaticBreakends(), somaticBreakendInterpreter::interpret))
                .driverSomaticBreakends(driverSomaticBreakends)
                .somaticHomozygousDisruptions(ConversionUtil.mapToIterable(linx.somaticHomozygousDisruptions(), LinxConversion::convert))
                .allGermlineStructuralVariants(allGermlineSvAnnotations)
                .otherGermlineBreakends(otherGermlineBreakends)
                .driverGermlineBreakends(driverGermlineBreakends)
                .driverSomaticDisruptions(driverSomaticDisruptions)
                .driverGermlineDisruptions(driverGermlineDisruptions)
                .germlineHomozygousDisruptions(ConversionUtil.mapToIterable(linx.germlineHomozygousDisruptions(), LinxConversion::convert))
                .build();
    }
}
