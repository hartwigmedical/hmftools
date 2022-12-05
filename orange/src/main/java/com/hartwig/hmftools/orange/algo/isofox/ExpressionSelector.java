package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.rna.GeneExpression;

import org.jetbrains.annotations.NotNull;

final class ExpressionSelector {

    private static final double HIGH_EXPRESSION_PERCENTILE_CUTOFF = 0.9;
    private static final double LOW_EXPRESSION_PERCENTILE_CUTOFF = 0.05;

    private ExpressionSelector() {
    }

    @NotNull
    public static List<GeneExpression> selectHighExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> oncogenes = extractGenesOfType(driverGenes, DriverCategory.ONCO);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            boolean isOncogene = oncogenes.contains(expression.geneName());
            boolean hasHighExpressionCohort = expression.percentileCohort() >= HIGH_EXPRESSION_PERCENTILE_CUTOFF;
            boolean hasHighExpressionCancer = expression.percentileCancer() >= HIGH_EXPRESSION_PERCENTILE_CUTOFF;
            if (isOncogene && hasHighExpressionCohort && hasHighExpressionCancer) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    public static List<GeneExpression> selectLowExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> suppressors = extractGenesOfType(driverGenes, DriverCategory.TSG);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            boolean isTSG = suppressors.contains(expression.geneName());
            boolean hasLowExpressionCohort = expression.percentileCohort() <= LOW_EXPRESSION_PERCENTILE_CUTOFF;
            boolean hasLowExpressionCancer = expression.percentileCancer() <= LOW_EXPRESSION_PERCENTILE_CUTOFF;
            if (isTSG && hasLowExpressionCohort && hasLowExpressionCancer) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    private static Set<String> extractGenesOfType(@NotNull List<DriverGene> driverGenes, @NotNull DriverCategory categoryToInclude) {
        Set<String> filtered = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.likelihoodType() == categoryToInclude) {
                filtered.add(driverGene.gene());
            }
        }
        return filtered;
    }
}
