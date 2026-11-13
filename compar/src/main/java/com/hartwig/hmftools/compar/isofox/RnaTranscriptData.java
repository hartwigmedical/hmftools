package com.hartwig.hmftools.compar.isofox;

import java.util.List;

import com.hartwig.hmftools.common.rna.TranscriptExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaTranscriptData extends ComparableItem
{
    public final TranscriptExpression Expression;

    public RnaTranscriptData(final TranscriptExpression expression, final List<FieldInfo> fields)
    {
        Expression = expression;

        addDoubleValue(RnaTranscriptDataComparer.Fields.AdjTPM.toString(), expression.tpm(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_TRANSCRIPT_DATA;
    }

    @Override
    public String key()
    {
        return Expression.transcriptName();
    }

    @Override
    public boolean reportable()
    {
        return false;
    }

    @Override
    public String geneName()
    {
        return Expression.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaTranscriptData otherData = (RnaTranscriptData)other;

        return otherData.Expression.transcriptName().equals(Expression.transcriptName());
    }
}
