package com.hartwig.hmftools.orange.report.pdfdata;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;

import org.jetbrains.annotations.Nullable;

public class SomaticFindingsDataFactory
{
    public static SomaticFindingsData build(
            final OrangeRecord report, @Nullable final OrangeConfig config, final PlotPathResolver plotPathResolver)
    {
        List<String> linxDriverPlotPaths = report.plots().linxDriverPlots().stream()
                .map(plotPathResolver::resolve)
                .collect(Collectors.toList());

        return new SomaticFindingsData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                report.tumorOnlyMode(),
                report.hasRna(),
                config != null && config.RnaSampleId != null,
                report.purple().somaticVariants(),
                report.purple().somaticGainsDels(),
                report.purple().armCopyNumberAbberations(),
                report.linx().fusions(),
                report.linx().somaticBreakends(),
                report.virusInterpreter(),
                report.sigAllocations(),
                linxDriverPlotPaths);
    }
}
