package com.hartwig.hmftools.compar.isofox;

import java.util.List;

import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaGeneData extends ComparableItem
{
    public final GeneExpression Expression;

    public RnaGeneData(final GeneExpression expression, final List<FieldInfo> fields)
    {
        Expression = expression;

        addIntValue(RnaGeneDataComparer.Fields.SplicedFragments.toString(), expression.splicedFragments(), fields);
        addIntValue(RnaGeneDataComparer.Fields.UnsplicedFragments.toString(), expression.unsplicedFragments(), fields);
        addDoubleValue(RnaGeneDataComparer.Fields.AdjTPM.toString(), expression.tpm(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_GENE_DATA;
    }

    @Override
    public String key()
    {
        return Expression.geneName();
    }

    @Override
    public boolean reportable()
    {
        return Expression.reportedStatus() == ReportedStatus.REPORTED;
    }

    @Override
    public String geneName()
    {
        return Expression.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaGeneData otherData = (RnaGeneData)other;

        return otherData.Expression.geneName().equals(Expression.geneName());
    }
}
