package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;

public class CuppaChapterDataFactory
{
    public static CuppaChapterData build(final OrangeRecord report, final PlotPathResolver plotPathResolver)
    {
        String cuppaSummaryPlotPath = report.plots().cuppaSummaryPlot() != null
                ? plotPathResolver.resolve(report.plots().cuppaSummaryPlot())
                : null;

        return new CuppaChapterData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                cuppaSummaryPlotPath);
    }
}
