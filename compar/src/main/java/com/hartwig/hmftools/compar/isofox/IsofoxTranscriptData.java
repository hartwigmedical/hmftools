package com.hartwig.hmftools.compar.isofox;

import com.hartwig.hmftools.common.rna.TranscriptExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record IsofoxTranscriptData(TranscriptExpression TranscriptExpression) implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_TRANSCRIPT_DATA;
    }

    @Override
    public String key()
    {
        return TranscriptExpression.transcriptName();
    }

    @Override
    public boolean reportable()
    {
        return false;
    }

    @Override
    public String geneName()
    {
        return TranscriptExpression.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final IsofoxTranscriptData otherData = (IsofoxTranscriptData)other;

        return otherData.TranscriptExpression.transcriptName().equals(TranscriptExpression.transcriptName());
    }
}
