package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.orange.algo.isofox.OrangeIsofoxTestFactory;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExpressionsTest
{
    @Test
    public void canCalcFoldChange()
    {
        GeneExpression normal = builder().tpm(1).medianTpmCancer(2D).medianTpmCohort(3D).build();
        assertEquals("0.5", Expressions.foldChangeType(normal));
        assertEquals("0.3", Expressions.foldChangeDatabase(normal));

        GeneExpression zeroMedian = builder().tpm(0).medianTpmCancer(0D).medianTpmCohort(0D).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeType(zeroMedian));
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeDatabase(zeroMedian));

        GeneExpression missingCancerTPM = builder().tpm(1).medianTpmCancer(null).medianTpmCohort(3D).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeType(missingCancerTPM));
    }

    @Test
    public void canFindExpressionByGene()
    {
        GeneExpression entry = builder().gene("gene 1").build();

        List<GeneExpression> expressions = Lists.newArrayList(entry);
        assertEquals(entry, Expressions.findByGene(expressions, "gene 1"));
        assertNull(Expressions.findByGene(expressions, "gene 2"));
    }

    @Test
    public void canFormatFieldsCorrectly()
    {
        GeneExpression expression = builder().tpm(1.234).percentileCancer(0.12345).percentileCohort(0.23445).build();

        assertEquals("1.2", Expressions.tpm(expression));
        assertEquals("0.12", Expressions.percentileType(expression));
        assertEquals("0.23", Expressions.percentileDatabase(expression));

        GeneExpression cancerExpressionMissing = builder().tpm(0).percentileCancer(null).percentileCohort(0).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.percentileType(cancerExpressionMissing));
    }

    @NotNull
    private static ImmutableGeneExpression.Builder builder()
    {
        return OrangeIsofoxTestFactory.geneExpressionBuilder();
    }
}