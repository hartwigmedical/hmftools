package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;

final class ExpressionSelector
{
    private static final double HIGH_EXPRESSION_PERCENTILE_CUTOFF = 0.9;
    private static final double LOW_EXPRESSION_PERCENTILE_CUTOFF = 0.05;

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
        return selectGenesMatchingCriteria(expressions,
                extractGenesOfType(driverGenes, DriverCategory.TSG),
                percentile -> percentile <= LOW_EXPRESSION_PERCENTILE_CUTOFF);
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
