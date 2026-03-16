package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ExpressionSelectorTest
{
    @Test
    public void testSelectReportableExpressionGenes()
    {
        com.hartwig.hmftools.common.rna.GeneExpression common = IsofoxTestFactory.geneExpressionBuilder().build();

        com.hartwig.hmftools.common.rna.GeneExpression highExpressionGene = create(GENE_NAME_1, 0.995, 0.95);
        com.hartwig.hmftools.common.rna.GeneExpression nonHighExpressionGene = create(GENE_NAME_2, 0.85, 0.95);
        com.hartwig.hmftools.common.rna.GeneExpression lowExpressionGene = create(GENE_NAME_3, 0.002, 0.95);

        List<com.hartwig.hmftools.common.rna.GeneExpression> geneGxpressions = Lists.newArrayList(
                highExpressionGene, nonHighExpressionGene, lowExpressionGene );

        Map<String,DriverGene> drivers = Maps.newHashMap();

        List<GeneExpression> highExpressionGenes = Lists.newArrayList();
        List<GeneExpression> lowExpressionGenes = Lists.newArrayList();

        IsofoxInterpreter.findExpressionOutliers(geneGxpressions, drivers, highExpressionGenes, lowExpressionGenes);

        assertTrue(highExpressionGenes.isEmpty());
        assertTrue(lowExpressionGenes.isEmpty());

        DriverGene driver1 = DriverGeneTestFactory.builder().gene(GENE_NAME_1).reportHighExpression(true).build();
        DriverGene driver2 = DriverGeneTestFactory.builder().gene(GENE_NAME_3).reportLowExpression(true).build();
        drivers.put(driver1.gene(), driver1);
        drivers.put(driver2.gene(), driver2);

        IsofoxInterpreter.findExpressionOutliers(geneGxpressions, drivers, highExpressionGenes, lowExpressionGenes);

        assertEquals(1, highExpressionGenes.size());
        assertTrue(highExpressionGenes.stream().anyMatch(x -> x.gene().equals(GENE_NAME_1)));

        assertEquals(1, lowExpressionGenes.size());
        assertTrue(lowExpressionGenes.stream().anyMatch(x -> x.gene().equals(GENE_NAME_3)));
    }

    private static com.hartwig.hmftools.common.rna.GeneExpression create(
            final String gene, double percentileCohort, @Nullable Double percentileCancer)
    {
        return IsofoxTestFactory.geneExpressionBuilder()
                .geneName(gene)
                .percentileCohort(percentileCohort)
                .percentileCancer(percentileCancer)
                .build();
    }
}