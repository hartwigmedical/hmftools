package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaNovelSpliceJunctionData extends ComparableItem
{
    public final String GeneName;
    public final String Chromosome;
    public final int JunctionStart;
    public final int JunctionEnd;

    public final AltSpliceJunctionType Type;
    public final AltSpliceJunctionContext RegionTypeStart;
    public final AltSpliceJunctionContext RegionTypeEnd;
    public final int FragmentCount;

    public RnaNovelSpliceJunctionData(
            final String geneName, final String chromosome, final int junctionStart, final int junctionEnd, final AltSpliceJunctionType type,
            final AltSpliceJunctionContext regionTypeStart, final AltSpliceJunctionContext regionTypeEnd, final int fragmentCount)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        JunctionStart = junctionStart;
        JunctionEnd = junctionEnd;
        Type = type;
        RegionTypeStart = regionTypeStart;
        RegionTypeEnd = regionTypeEnd;
        FragmentCount = fragmentCount;
    }

    public static RnaNovelSpliceJunctionData from(final NovelSpliceJunction novelSpliceJunction, final List<FieldInfo> fields)
    {
        RnaNovelSpliceJunctionData junction = new RnaNovelSpliceJunctionData(
                novelSpliceJunction.geneName(), novelSpliceJunction.chromosome(), novelSpliceJunction.junctionStart(),
                novelSpliceJunction.junctionEnd(), novelSpliceJunction.type(),
                novelSpliceJunction.regionStart(), novelSpliceJunction.regionEnd(),
                novelSpliceJunction.fragmentCount());

        junction.addAllValues(fields);

        return junction;
    }

    private void addAllValues(final List<FieldInfo> fields)
    {
        addStringValue(RnaNovelSpliceJunctionComparer.Fields.Type.toString(), Type.toString(), fields);
        addIntValue(RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString(), FragmentCount, fields);
        addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString(), RegionTypeStart.toString(), fields);
        addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString(), RegionTypeEnd.toString(), fields);
    }

    public static RnaNovelSpliceJunctionData fromTruthset(final List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        String key = truthsetValues.get(0).Key;
        String[] keyParts = key.split(":", 4);

        String gene = keyParts[0];
        String chromosome = keyParts[1];
        int positionUp = Integer.parseInt(keyParts[2]);
        int positionDown = Integer.parseInt(keyParts[3]);

        AltSpliceJunctionType type = AltSpliceJunctionType.UNKNOWN;
        AltSpliceJunctionContext junctionTypeUp = AltSpliceJunctionContext.UNKNOWN;
        AltSpliceJunctionContext junctionTypeDown = AltSpliceJunctionContext.UNKNOWN;
        int fragmentCount = 0;

        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.Type.toString()))
                type = AltSpliceJunctionType.valueOf(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString()))
                fragmentCount = Integer.parseInt(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString()))
                junctionTypeUp = AltSpliceJunctionContext.valueOf(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString()))
                junctionTypeDown = AltSpliceJunctionContext.valueOf(truthsetValue.Value);
        }

        RnaNovelSpliceJunctionData junction = new RnaNovelSpliceJunctionData(
                gene, chromosome, positionUp, positionDown, type,
                junctionTypeUp, junctionTypeDown, fragmentCount);

        junction.addTruthsetValues(truthsetValues, fields);
        return junction;
    }

    private void addTruthsetValues(List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.Type.toString()))
                addStringValue(RnaNovelSpliceJunctionComparer.Fields.Type.toString(), Type.toString(), fields);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString()))
                addIntValue(RnaNovelSpliceJunctionComparer.Fields.FragmentCount.toString(), FragmentCount, fields);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString()))
                addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionStart.toString(), RegionTypeStart.toString(), fields);
            else if(truthsetValue.FieldName.equals(RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString()))
                addStringValue(RnaNovelSpliceJunctionComparer.Fields.RegionEnd.toString(), RegionTypeEnd.toString(), fields);
        }
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_NOVEL_SPLICE_JUNCTION;
    }

    @Override
    public String key()
    {
        return format("%s %s:%d-%d", GeneName, Chromosome, JunctionStart, JunctionEnd);
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public String geneName()
    {
        return GeneName;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaNovelSpliceJunctionData otherData = (RnaNovelSpliceJunctionData)other;

        if(!otherData.GeneName.equals(GeneName))
            return false;

        if(!otherData.Chromosome.equals(Chromosome))
            return false;

        if(otherData.JunctionStart != JunctionStart)
            return false;

        return otherData.JunctionEnd == JunctionEnd;
    }
}
