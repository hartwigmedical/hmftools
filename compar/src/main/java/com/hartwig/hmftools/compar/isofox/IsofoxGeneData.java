package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record IsofoxGeneData(GeneExpression GeneExpression) implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_GENE_DATA;
    }

    @Override
    public String key()
    {
        return GeneExpression.geneName();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%d", GeneExpression.splicedFragments()));
        values.add(format("%d", GeneExpression.unsplicedFragments()));
        values.add(format("%.2f", GeneExpression.tpm()));
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
        return GeneExpression.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final IsofoxGeneData otherData = (IsofoxGeneData)other;

        return otherData.GeneExpression.geneName().equals(GeneExpression.geneName());
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final GeneExpression ref = GeneExpression;
        final GeneExpression otherData = ((IsofoxGeneData) other).GeneExpression;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_SPLICED_FRAGS, ref.splicedFragments(), otherData.splicedFragments(), thresholds);
        checkDiff(diffs, FLD_UNSPLICED_FRAGS, ref.unsplicedFragments(), otherData.unsplicedFragments(), thresholds);
        checkDiff(diffs, FLD_ADJ_TPM, ref.tpm(), otherData.tpm(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
