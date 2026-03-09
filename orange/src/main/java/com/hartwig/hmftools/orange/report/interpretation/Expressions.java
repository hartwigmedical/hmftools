package com.hartwig.hmftools.orange.report.interpretation;

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
}
