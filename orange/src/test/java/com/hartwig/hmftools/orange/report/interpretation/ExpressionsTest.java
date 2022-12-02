package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExpressionsTest {

    @Test
    public void canCalcFoldChange() {
        GeneExpression normal = builder().tpm(1).medianTpmCancer(2).medianTpmCohort(3).build();
        assertEquals("0.5", Expressions.foldChangeType(normal));
        assertEquals("0.3", Expressions.foldChangeDatabase(normal));

        GeneExpression zeroMedian = builder().tpm(0).medianTpmCancer(0).medianTpmCohort(0).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeType(zeroMedian));
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.foldChangeDatabase(zeroMedian));
    }

    @Test
    public void canFindExpressionByGene() {
        GeneExpression entry = builder().geneName("gene 1").build();

        List<GeneExpression> expressions = Lists.newArrayList(entry);
        assertEquals(entry, Expressions.findByGene(expressions, "gene 1"));
        assertNull(Expressions.findByGene(expressions, "gene 2"));
    }

    @Test
    public void canFormatFieldsCorrectly() {
        GeneExpression expression = builder().tpm(1.234).percentileCancer(0.12345).percentileCohort(0.23445).build();

        assertEquals("1.2", Expressions.tpm(expression));
        assertEquals("0.12", Expressions.percentileType(expression));
        assertEquals("0.23", Expressions.percentileDatabase(expression));
    }

    @NotNull
    private static ImmutableGeneExpression.Builder builder() {
        return IsofoxTestFactory.geneExpressionBuilder();
    }
}