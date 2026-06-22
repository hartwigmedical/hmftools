package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;

public class RnaFindingsDataFactory
{
    public static RnaFindingsData build(final IsofoxRecord isofox)
    {
        return new RnaFindingsData(
                QcStatusInterpretation.hasRnaFail(isofox),
                isofox.summary(),
                isofox.highExpressionGenes(),
                isofox.fusions(),
                isofox.novelSpliceJunctions());
    }
}
