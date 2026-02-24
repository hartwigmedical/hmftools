package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.report.ReportResources.formatFoldChangeField;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.jetbrains.annotations.Nullable;

public final class Expressions
{
    @Nullable
    public static GeneExpression findByGene(final List<GeneExpression> expressions, final String geneToFind)
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

    public static String formatFoldChange(final GeneExpression expression)
    {
        return toFoldChangeDisplay(expression.tpm(), expression.medianTpmCohort());
    }

    public static String formatFoldChangeCancer(final GeneExpression expression)
    {
        return toFoldChangeDisplay(expression.tpm(), expression.medianTpmCancer());
    }

    private static String toFoldChangeDisplay(double expression, @Nullable Double median)
    {
        if(median == null || Doubles.lessOrEqual(median, 0))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        double foldChange = expression / median;
        return formatFoldChangeField(foldChange);
    }


    /*
    public static String tpm(final GeneExpression expression)
    {
        return formatSingleDigitDecimal(expression.tpm());
    }

    public static String percentileType(final GeneExpression expression)
    {
        return expression.percentileCancer() == null ? ReportResources.NOT_AVAILABLE : formatTwoDigitDecimal(expression.percentileCancer());
    }

    public static String percentileDatabase(final GeneExpression expression)
    {
        return formatTwoDigitDecimal(expression.percentileCohort());
    }

    public static String foldChangeDatabase(final GeneExpression expression)
    {
        return toFoldChange(expression.tpm(), expression.medianTpmCohort());
    }

    private static String toFoldChange(double expression, @Nullable Double median)
    {
        if(median == null || Doubles.lessOrEqual(median, 0))
        {
            return ReportResources.NOT_AVAILABLE;
        }

        double foldChange = expression / median;
        return foldChange > 1000 ? ">1000" : formatSingleDigitDecimal(foldChange);
    }
    */
}
