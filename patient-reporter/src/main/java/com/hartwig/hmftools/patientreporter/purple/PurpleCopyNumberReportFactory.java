package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReportType;
import com.hartwig.hmftools.patientreporter.copynumber.ImmutableCopyNumberReport;

import org.jetbrains.annotations.NotNull;

class PurpleCopyNumberReportFactory {

    @VisibleForTesting
    static final double ABS_GAIN = 8;
    @VisibleForTesting
    static final double REL_GAIN = 2.2;
    @VisibleForTesting
    static final double ABS_LOSS = 0.5;

    @NotNull
    static List<CopyNumberReport> createReport(final double samplePloidy, @NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        final List<CopyNumberReport> result = Lists.newArrayList();
        for (final GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            final CopyNumberReportType type = type(samplePloidy, geneCopyNumber.minCopyNumber());
            if (type != CopyNumberReportType.NEUTRAL) {
                result.add(ImmutableCopyNumberReport.builder()
                        .chromosome(geneCopyNumber.chromosome())
                        .chromosomeBand(geneCopyNumber.chromosomeBand())
                        .gene(geneCopyNumber.gene())
                        .copyNumber(geneCopyNumber.value())
                        .type(type)
                        .build());
            }
        }
        return result;
    }

    @VisibleForTesting
    @NotNull
    static CopyNumberReportType type(final double samplePloidy, final double copyNumber) {
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            return CopyNumberReportType.LOSS;
        }

        if (Doubles.greaterOrEqual(copyNumber, ABS_GAIN)) {
            return CopyNumberReportType.GAIN;
        }

        double relativeCopyNumber = copyNumber / samplePloidy;
        if (Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN)) {
            return CopyNumberReportType.GAIN;
        }

        return CopyNumberReportType.NEUTRAL;
    }
}
