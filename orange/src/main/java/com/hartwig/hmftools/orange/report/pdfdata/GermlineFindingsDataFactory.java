package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;

public class GermlineFindingsDataFactory
{
    public static GermlineFindingsData build(final OrangeRecord report)
    {
        return new GermlineFindingsData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                report.referenceId() != null,
                report.hasRna(),
                report.purple().germlineDrivers(),
                report.purple().germlineVariants(),
                report.purple().germlineGainsDels(),
                report.linx().germlineBreakends(),
                report.purple().fit().qc().germlineAberrations(),
                report.peach());
    }
}
