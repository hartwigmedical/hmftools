package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.common.rna.GeneExpression;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExpressionSelectorTest {

    @Test
    public void canSelectHighExpressionGenes() {
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
    public void canSelectLowExpressionGenes() {
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

    @NotNull
    private static GeneExpression create(@NotNull String gene, double percentileCohort, double percentileCancer) {
        return IsofoxTestFactory.geneExpressionBuilder()
                .geneName(gene)
                .percentileCohort(percentileCohort)
                .percentileCancer(percentileCancer)
                .build();
    }
}