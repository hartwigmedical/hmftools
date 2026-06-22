package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;

public class QualityControlDataFactory
{
    public static QualityControlData build(final OrangeRecord report, final PlotPathResolver plotPathResolver)
    {
        String qSeePlotPath = report.plots().qSeePlot() != null
                ? plotPathResolver.resolve(report.plots().qSeePlot())
                : null;

        return new QualityControlData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                qSeePlotPath);
    }
}
