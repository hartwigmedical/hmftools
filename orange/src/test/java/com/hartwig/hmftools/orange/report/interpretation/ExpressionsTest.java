package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;

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
        assertEquals("0.3", Expressions.formatFoldChange(normal)); // uses pan-cancer values for now

        GeneExpression zeroMedian = builder().tpm(0).medianTpmCancer(0D).medianTpmCohort(0D).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.formatFoldChangeCancer(zeroMedian));

        GeneExpression missingCancerTPM = builder().tpm(1).medianTpmCancer(null).medianTpmCohort(3D).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Expressions.formatFoldChangeCancer(missingCancerTPM));
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

        assertEquals("1.2", formatTpmField(expression.tpm()));
        assertEquals("0.12", formatPercentileField(expression.percentileCancer()));
        assertEquals("0.23", formatPercentileField(expression.percentileCohort()));
    }

    @NotNull
    private static ImmutableGeneExpression.Builder builder()
    {
        return OrangeIsofoxTestFactory.geneExpressionBuilder();
    }
}