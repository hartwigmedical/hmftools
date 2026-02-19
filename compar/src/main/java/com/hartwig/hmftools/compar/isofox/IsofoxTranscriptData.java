package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.TranscriptExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record IsofoxTranscriptData(TranscriptExpression TranscriptExpression) implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_TRANSCRIPT_DATA;
    }

    @Override
    public String key()
    {
        return TranscriptExpression.transcriptName();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", TranscriptExpression.geneName()));
        values.add(format("%.2f", TranscriptExpression.tpm()));
        return values;
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

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final TranscriptExpression ref = TranscriptExpression;
        final TranscriptExpression otherData = ((IsofoxTranscriptData) other).TranscriptExpression;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_GENE_NAME, ref.geneName(), otherData.geneName());
        checkDiff(diffs, FLD_ADJ_TPM, ref.tpm(), otherData.tpm(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
