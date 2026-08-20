package com.hartwig.hmftools.compar.isofox;

import java.util.List;

import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaNovelSpliceJunctionData extends ComparableItem
{
    public final NovelSpliceJunction Junction;

    public RnaNovelSpliceJunctionData(final NovelSpliceJunction novelSpliceJunction, final List<FieldInfo> fields)
    {
        Junction = novelSpliceJunction;

        addStringValue(RnaNovelSpliceJunctionComparer.Fields.Type.toString(), novelSpliceJunction.type().toString(), fields);
        addIntValue(RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString(), novelSpliceJunction.fragmentCount(), fields);
        addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString(), novelSpliceJunction.regionStart().toString(), fields);
        addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString(), novelSpliceJunction.regionEnd().toString(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_NOVEL_SPLICE_JUNCTION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%d", Junction.geneName(), Junction.chromosome(),
                Junction.junctionStart(), Junction.junctionEnd());

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
        final RnaNovelSpliceJunctionData otherData = (RnaNovelSpliceJunctionData)other;

        if(!otherData.Junction.geneName().equals(Junction.geneName())){
            return false;
        }

        if(!otherData.Junction.chromosome().equals(Junction.chromosome()))
            return false;

        if(otherData.Junction.junctionStart() != Junction.junctionStart())
            return false;

        return otherData.Junction.junctionEnd() == Junction.junctionEnd();
    }
}
