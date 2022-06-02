package com.hartwig.hmftools.orange.report.interpretation;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Expressions {

    private static final Logger LOGGER = LogManager.getLogger(Expressions.class);

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat TWO_DIGIT = ReportResources.decimalFormat("#.##");

    private Expressions() {
    }

    @Nullable
    public static GeneExpression findByGene(@NotNull List<GeneExpression> expressions, @NotNull String geneToFind) {
        for (GeneExpression expression : expressions) {
            if (expression.geneName().equals(geneToFind)) {
                return expression;
            }
        }

        LOGGER.warn("Could not find expression data for gene '{}'", geneToFind);
        return null;
    }

    @NotNull
    public static String tpm(@NotNull GeneExpression expression) {
        return SINGLE_DIGIT.format(expression.tpm());
    }

    @NotNull
    public static String percentileType(@NotNull GeneExpression expression) {
        return TWO_DIGIT.format(expression.percentileCancer());
    }

    @NotNull
    public static String foldChangeType(@NotNull GeneExpression expression) {
        return toFoldChange(expression.tpm(), expression.medianTpmCancer());
    }

    @NotNull
    public static String percentileDatabase(@NotNull GeneExpression expression) {
        return TWO_DIGIT.format(expression.percentileCohort());
    }

    @NotNull
    public static String foldChangeDatabase(@NotNull GeneExpression expression) {
        return toFoldChange(expression.tpm(), expression.medianTpmCohort());
    }

    @NotNull
    private static String toFoldChange(double expression, double median) {
        if (Doubles.lessOrEqual(median, 0)) {
            return ReportResources.NOT_AVAILABLE;
        }

        double foldChange = expression / median;
        return foldChange > 1000 ? ">1000" : SINGLE_DIGIT.format(foldChange);
    }
}
