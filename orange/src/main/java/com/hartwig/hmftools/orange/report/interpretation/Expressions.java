package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Expressions
{
    @Nullable
    public static GeneExpression findByGene(@NotNull List<GeneExpression> expressions, @NotNull String geneToFind)
    {
        for(GeneExpression expression : expressions)
        {
            if(expression.gene().equals(geneToFind))
            {
                return expression;
            }
        }

        LOGGER.warn("Could not find expression data for gene '{}'", geneToFind);
        return null;
    }

    @NotNull
    public static String tpm(@NotNull GeneExpression expression)
    {
        return formatSingleDigitDecimal(expression.tpm());
    }

    @NotNull
    public static String percentileType(@NotNull GeneExpression expression)
    {
        return expression.percentileCancer() == null ? ReportResources.NOT_AVAILABLE : formatTwoDigitDecimal(expression.percentileCancer());
    }

    @NotNull
    public static String foldChangeType(@NotNull GeneExpression expression)
    {
        return toFoldChange(expression.tpm(), expression.medianTpmCancer());
    }

    @NotNull
    public static String percentileDatabase(@NotNull GeneExpression expression)
    {
        return formatTwoDigitDecimal(expression.percentileCohort());
    }

    @NotNull
    public static String foldChangeDatabase(@NotNull GeneExpression expression)
    {
        return toFoldChange(expression.tpm(), expression.medianTpmCohort());
    }

    @NotNull
    private static String toFoldChange(double expression, @Nullable Double median)
    {
        if(median == null || Doubles.lessOrEqual(median, 0))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        double foldChange = expression / median;
        return foldChange > 1000 ? ">1000" : formatSingleDigitDecimal(foldChange);
    }
}
