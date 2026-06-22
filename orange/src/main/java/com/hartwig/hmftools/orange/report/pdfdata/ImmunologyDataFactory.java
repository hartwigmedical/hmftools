package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;

public class ImmunologyDataFactory
{
    public static ImmunologyData build(final OrangeRecord report)
    {
        return new ImmunologyData(
                QcStatusInterpretation.hasPurpleFail(report.purple().fit().qc()),
                report.lilac(),
                report.hasRna(),
                report.experimentType() == ExperimentType.WHOLE_GENOME);
    }
}
