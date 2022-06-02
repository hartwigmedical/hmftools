package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.junit.Test;

public class ExpressionsTest {

    @Test
    public void canCalcFoldChange() {
        GeneExpression normal = IsofoxTestFactory.geneExpressionBuilder().tpm(1).medianTpmCancer(2).build();
        assertEquals("0.5", Expressions.foldChangeType(normal));

        GeneExpression zeroMedian = IsofoxTestFactory.geneExpressionBuilder().tpm(0).medianTpmCancer(0).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeType(zeroMedian));
    }
}