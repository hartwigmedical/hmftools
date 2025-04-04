package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ExpressionSelectorTest
{
    @Test
    public void shouldSelectHighExpressionGenes()
    {
        GeneExpression highExpressionGene1 = create("gene 1", 0.95, 0.95);
        GeneExpression nonHighExpressionGene1 = create("gene 1", 0.85, 0.95);
        GeneExpression highExpressionGene2 = create("gene 2", 0.95, 0.95);
        List<GeneExpression> expressions = Lists.newArrayList(highExpressionGene1, nonHighExpressionGene1, highExpressionGene2);

        DriverGene driver1 = DriverGeneTestFactory.builder().gene("gene 1").likelihoodType(DriverCategory.ONCO).build();
        DriverGene driver2 = DriverGeneTestFactory.builder().gene("gene 2").likelihoodType(DriverCategory.TSG).build();
        List<DriverGene> drivers = Lists.newArrayList(driver1, driver2);

        List<GeneExpression> highExpression = ExpressionSelector.selectHighExpressionGenes(expressions, drivers);

        assertEquals(1, highExpression.size());
        assertTrue(highExpression.contains(highExpressionGene1));
    }

    @Test
    public void shouldSelectHighExpressionGeneWithoutCohort()
    {
        GeneExpression highExpressionGene = create("gene 1", 0.95, null);
        DriverGene driver = DriverGeneTestFactory.builder().gene("gene 1").likelihoodType(DriverCategory.ONCO).build();

        List<GeneExpression> highExpression = ExpressionSelector.selectHighExpressionGenes(List.of(highExpressionGene), List.of(driver));

        assertEquals(1, highExpression.size());
        assertTrue(highExpression.contains(highExpressionGene));
    }

    @Test
    public void shouldSelectLowExpressionGenes()
    {
        GeneExpression lowExpressionGene1 = create("gene 1", 0.02, 0.02);
        GeneExpression nonLowExpressionGene1 = create("gene 1", 0.02, 0.08);
        GeneExpression lowExpressionGene2 = create("gene 2", 0.02, 0.02);
        List<GeneExpression> expressions = Lists.newArrayList(lowExpressionGene1, nonLowExpressionGene1, lowExpressionGene2);

        DriverGene driver1 = DriverGeneTestFactory.builder().gene("gene 1").likelihoodType(DriverCategory.ONCO).build();
        DriverGene driver2 = DriverGeneTestFactory.builder().gene("gene 2").likelihoodType(DriverCategory.TSG).build();
        List<DriverGene> drivers = Lists.newArrayList(driver1, driver2);

        List<GeneExpression> lowExpression = ExpressionSelector.selectLowExpressionGenes(expressions, drivers);

        assertEquals(1, lowExpression.size());
        assertTrue(lowExpression.contains(lowExpressionGene2));
    }

    @Test
    public void shouldSelectLowExpressionGenesWithoutCohort()
    {
        GeneExpression lowExpressionGene = create("gene 2", 0.02, null);
        DriverGene driver = DriverGeneTestFactory.builder().gene("gene 2").likelihoodType(DriverCategory.TSG).build();

        List<GeneExpression> lowExpression = ExpressionSelector.selectLowExpressionGenes(List.of(lowExpressionGene), List.of(driver));

        assertEquals(1, lowExpression.size());
        assertTrue(lowExpression.contains(lowExpressionGene));
    }

    @NotNull
    private static GeneExpression create(@NotNull String gene, double percentileCohort, @Nullable Double percentileCancer)
    {
        return OrangeIsofoxTestFactory.geneExpressionBuilder()
                .gene(gene)
                .percentileCohort(percentileCohort)
                .percentileCancer(percentileCancer)
                .build();
    }
}