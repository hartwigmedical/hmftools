package com.hartwig.hmftools.compar.isofox;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

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
}
