package com.hartwig.hmftools.compar.isofox;

import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record IsofoxGeneData(GeneExpression GeneExpression) implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_GENE_DATA;
    }

    @Override
    public String key()
    {
        return GeneExpression.geneName();
    }

    @Override
    public boolean reportable()
    {
        return GeneExpression.reportedStatus() == ReportedStatus.REPORTED;
    }

    @Override
    public String geneName()
    {
        return GeneExpression.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final IsofoxGeneData otherData = (IsofoxGeneData)other;

        return otherData.GeneExpression.geneName().equals(GeneExpression.geneName());
    }
}
