package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.HIGH_EXPRESSION_PERCENTILE_CUTOFF;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.LOW_EXPRESSION_PERCENTILE_CUTOFF;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.MAX_EXPRESSION_GENE_COUNT;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;

final class ExpressionSelector
{
    public static List<GeneExpression> selectHighExpressionGenes(
            final List<GeneExpression> expressions, final List<DriverGene> driverGenes)
    {
        return selectGenesMatchingCriteria(expressions,
                extractGenesOfType(driverGenes, DriverCategory.ONCO),
                percentile -> percentile >= HIGH_EXPRESSION_PERCENTILE_CUTOFF);
    }

    public static List<GeneExpression> selectLowExpressionGenes(
            final List<GeneExpression> expressions, final List<DriverGene> driverGenes)
    {
        List<GeneExpression> lowGeneExpressions = selectGenesMatchingCriteria(expressions,
                extractGenesOfType(driverGenes, DriverCategory.TSG),
                percentile -> percentile <= LOW_EXPRESSION_PERCENTILE_CUTOFF);

        if(lowGeneExpressions.size() > MAX_EXPRESSION_GENE_COUNT)
            lowGeneExpressions = lowGeneExpressions.subList(0, MAX_EXPRESSION_GENE_COUNT);

        return lowGeneExpressions;
    }

    private static List<GeneExpression> selectGenesMatchingCriteria(
            final List<GeneExpression> expressions, Set<String> genesOfType, Predicate<Double> evaluatePercentileThreshold)
    {
        return expressions.stream().filter(expression ->
        {
            boolean geneMatchesType = genesOfType.contains(expression.gene());
            boolean percentileCancerMeetsThreshold =
                    expression.percentileCancer() == null || evaluatePercentileThreshold.test(expression.percentileCancer());
            boolean percentileCohortMeetsThreshold = evaluatePercentileThreshold.test(expression.percentileCohort());
            return geneMatchesType && percentileCancerMeetsThreshold && percentileCohortMeetsThreshold;
        }).collect(Collectors.toList());
    }

    private static Set<String> extractGenesOfType(final List<DriverGene> driverGenes, final DriverCategory categoryToInclude)
    {
        return driverGenes.stream()
                .filter(driverGene -> driverGene.likelihoodType() == categoryToInclude)
                .map(DriverGene::gene)
                .collect(Collectors.toSet());
    }
}
