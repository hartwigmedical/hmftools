package com.hartwig.hmftools.orange.report.interpretation;

import java.text.DecimalFormat;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.jetbrains.annotations.NotNull;

public final class Expressions {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private Expressions() {
    }

    @NotNull
    public static String tpm(@NotNull GeneExpression expression) {
        return SINGLE_DIGIT.format(expression.tpm());
    }

    @NotNull
    public static String percentileType(@NotNull GeneExpression expression) {
        return PERCENTAGE.format(expression.percentileCancer() * 100);
    }

    @NotNull
    public static String foldChangeType(@NotNull GeneExpression expression) {
        return formatFoldChange(expression.tpm() / expression.medianTpmCancer());
    }

    @NotNull
    public static String percentileDatabase(@NotNull GeneExpression expression) {
        return PERCENTAGE.format(expression.percentileCohort() * 100);
    }

    @NotNull
    public static String foldChangeDatabase(@NotNull GeneExpression expression) {
        return formatFoldChange(expression.tpm() / expression.medianTpmCohort());
    }

    @NotNull
    private static String formatFoldChange(double foldChange) {
        return foldChange > 1000 ? ">1000" : SINGLE_DIGIT.format(foldChange);
    }
}
