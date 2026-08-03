package com.hartwig.hmftools.compar.isofox;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class NovelSpliceJunctionData extends ComparableItem
{
    public final NovelSpliceJunction Junction;
    public final BasePosition ComparisonPositionStart;
    public final BasePosition ComparisonPositionEnd;

    public NovelSpliceJunctionData(
            final NovelSpliceJunction novelSpliceJunction, final BasePosition comparisonPositionStart,
            final BasePosition comparisonPositionEnd, final List<FieldInfo> fields)
    {
        Junction = novelSpliceJunction;
        ComparisonPositionStart = comparisonPositionStart;
        ComparisonPositionEnd = comparisonPositionEnd;

        addStringValue(NovelSpliceJunctionComparer.Fields.Type.toString(), novelSpliceJunction.type().toString(), fields);
        addIntValue(NovelSpliceJunctionComparer.Fields.FragmentCount.toString(), novelSpliceJunction.fragmentCount(), fields);
        addStringValue(NovelSpliceJunctionComparer.Fields.RegionStart.toString(), novelSpliceJunction.regionStart().toString(), fields);
        addStringValue(NovelSpliceJunctionComparer.Fields.RegionEnd.toString(), novelSpliceJunction.regionEnd().toString(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.NOVEL_SPLICE_JUNCTION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%d", Junction.geneName(), Junction.chromosome(),
                Junction.junctionStart(), Junction.junctionEnd());

        boolean startLifted = ComparisonPositionStart.Position != Junction.junctionStart()
                || !ComparisonPositionStart.Chromosome.equals(Junction.chromosome());

        boolean endLifted = ComparisonPositionEnd.Position != Junction.junctionEnd()
                || !ComparisonPositionEnd.Chromosome.equals(Junction.chromosome());

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
        return Junction.geneName();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final NovelSpliceJunctionData otherData = (NovelSpliceJunctionData)other;

        if(!otherData.Junction.geneName().equals(Junction.geneName())){
            return false;
        }
        if(!otherData.Junction.chromosome().equals(ComparisonPositionStart.Chromosome))
        {
            return false;
        }
        if(!otherData.Junction.chromosome().equals(ComparisonPositionEnd.Chromosome))
        {
            return false;
        }
        if(otherData.Junction.junctionStart() != ComparisonPositionStart.Position)
        {
            return false;
        }
        return otherData.Junction.junctionEnd() == ComparisonPositionEnd.Position;
    }
}
