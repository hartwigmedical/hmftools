package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record NovelSpliceJunctionData(NovelSpliceJunction NovelSpliceJunction, BasePosition ComparisonPositionStart,
                                      BasePosition ComparisonPositionEnd) implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.NOVEL_SPLICE_JUNCTION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%d", NovelSpliceJunction.geneName(), NovelSpliceJunction.chromosome(),
                NovelSpliceJunction.junctionStart(), NovelSpliceJunction.junctionEnd());

        boolean startLifted = ComparisonPositionStart.Position != NovelSpliceJunction.junctionStart()
                || !ComparisonPositionStart.Chromosome.equals(NovelSpliceJunction.chromosome());

        boolean endLifted = ComparisonPositionEnd.Position != NovelSpliceJunction.junctionEnd()
                || !ComparisonPositionEnd.Chromosome.equals(NovelSpliceJunction.chromosome());

        if(startLifted || endLifted)
            key += String.format(" liftover(%s-%s)", ComparisonPositionStart, ComparisonPositionEnd);

        return key;
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", NovelSpliceJunction.type()));
        values.add(format("%d", NovelSpliceJunction.fragmentCount()));
        values.add(format("%s", NovelSpliceJunction.regionStart()));
        values.add(format("%s", NovelSpliceJunction.regionEnd()));

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
        return NovelSpliceJunction.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final NovelSpliceJunctionData otherData = (NovelSpliceJunctionData)other;

        if(!otherData.NovelSpliceJunction.geneName().equals(NovelSpliceJunction.geneName())){
            return false;
        }
        if(!otherData.NovelSpliceJunction.chromosome().equals(ComparisonPositionStart.Chromosome))
        {
            return false;
        }
        if(!otherData.NovelSpliceJunction.chromosome().equals(ComparisonPositionEnd.Chromosome))
        {
            return false;
        }
        if(otherData.NovelSpliceJunction.junctionStart() != ComparisonPositionStart.Position)
        {
            return false;
        }
        return otherData.NovelSpliceJunction.junctionEnd() == ComparisonPositionEnd.Position;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final NovelSpliceJunction ref = NovelSpliceJunction;
        final NovelSpliceJunction otherData = ((NovelSpliceJunctionData) other).NovelSpliceJunction;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_ALT_SJ_TYPE, ref.type().toString(), otherData.type().toString());
        checkDiff(diffs, FLD_FRAG_COUNT, ref.fragmentCount(), otherData.fragmentCount(), thresholds);
        checkDiff(diffs, FLD_REGION_START, ref.regionStart().toString(), otherData.regionStart().toString());
        checkDiff(diffs, FLD_REGION_END, ref.regionEnd().toString(), otherData.regionEnd().toString());

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
