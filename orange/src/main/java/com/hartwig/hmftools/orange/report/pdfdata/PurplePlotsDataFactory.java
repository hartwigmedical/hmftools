package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;

public class PurplePlotsDataFactory
{
    public static PurplePlotsData build(final OrangeRecord report, final PlotPathResolver plotPathResolver)
    {
        return new PurplePlotsData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                plotPathResolver.resolve(report.plots().purpleInputCircosPlot()),
                plotPathResolver.resolve(report.plots().purpleCopyNumberPlot()),
                plotPathResolver.resolve(report.plots().purpleClonalityPlot()),
                plotPathResolver.resolve(report.plots().purplePurityRangePlot()),
                plotPathResolver.resolve(report.plots().purpleMinorAlleleMapPlot()),
                plotPathResolver.resolve(report.plots().purpleRainfallPlot()));
    }
}
