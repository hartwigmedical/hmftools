package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.formFusionName;

import java.util.List;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaFusionData extends ComparableItem
{
    public final String Name;
    public final String ChromosomeUp;
    public final String ChromosomeDown;
    public final int PositionUp;
    public final int PositionDown;
    public final KnownFusionType KnownType;
    public final String JunctionTypeUp;
    public final String JunctionTypeDown;
    public final int SplitFragments;

    public RnaFusionData(
            final String name, final String chromosomeUp, final String chromosomeDown, final int positionUp,
            final int positionDown, final KnownFusionType knownType,
            final String junctionTypeUp, final String junctionTypeDown, final int splitFragments)
    {
        Name = name;
        ChromosomeUp = chromosomeUp;
        ChromosomeDown = chromosomeDown;
        PositionUp = positionUp;
        PositionDown = positionDown;
        KnownType = knownType;
        JunctionTypeUp = junctionTypeUp;
        JunctionTypeDown = junctionTypeDown;
        SplitFragments = splitFragments;
    }

    public static RnaFusionData from(final RnaFusion rnaFusion, final List<FieldInfo> fields)
    {
        RnaFusionData fusionData = new RnaFusionData(
                rnaFusion.name(), rnaFusion.chromosomeUp(), rnaFusion.chromosomeDown(), rnaFusion.positionUp(), rnaFusion.positionDown(),
                rnaFusion.knownType(), rnaFusion.junctionTypeUp(), rnaFusion.junctionTypeDown(), rnaFusion.splitFragments());

        fusionData.addAllValues(fields);
        return fusionData;
    }

    private void addAllValues(final List<FieldInfo> fields)
    {
        addStringValue(RnaFusionComparer.Fields.KnownType.toString(), KnownType.toString(), fields);
        addStringValue(RnaFusionComparer.Fields.JuncTypeUp.toString(), JunctionTypeUp, fields);
        addStringValue(RnaFusionComparer.Fields.JuncTypeDown.toString(), JunctionTypeDown, fields);
        addIntValue(RnaFusionComparer.Fields.SplitFrags.toString(), SplitFragments, fields);
    }

    public static RnaFusionData fromTruthset(final List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        String key = truthsetValues.get(0).Key;
        String[] keyParts = key.split(":", 6);

        String geneUp = keyParts[0];
        String geneDown = keyParts[1];
        String name = formFusionName(geneUp, geneDown);
        String chromosomeUp = keyParts[2];
        int positionUp = Integer.parseInt(keyParts[3]);
        String chromosomeDown = keyParts[4];
        int positionDown = Integer.parseInt(keyParts[5]);

        KnownFusionType knownType = KnownFusionType.NONE;
        String junctionTypeUp = "";
        String junctionTypeDown = "";
        int splitFragments = 0;

        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.KnownType.toString()))
                knownType = KnownFusionType.valueOf(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.SplitFrags.toString()))
                splitFragments = Integer.parseInt(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.JuncTypeUp.toString()))
                junctionTypeUp = truthsetValue.Value;
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.JuncTypeDown.toString()))
                junctionTypeDown = truthsetValue.Value;
        }

        RnaFusionData fusionData = new RnaFusionData(
                name, chromosomeUp, chromosomeDown, positionUp, positionDown, knownType, junctionTypeUp, junctionTypeDown, splitFragments);

        fusionData.addTruthsetValues(truthsetValues, fields);
        return fusionData;
    }

    private void addTruthsetValues(List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.KnownType.toString()))
                addStringValue(RnaFusionComparer.Fields.KnownType.toString(), KnownType.toString(), fields);
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.SplitFrags.toString()))
                addIntValue(RnaFusionComparer.Fields.SplitFrags.toString(), SplitFragments, fields);
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.JuncTypeUp.toString()))
                addStringValue(RnaFusionComparer.Fields.JuncTypeUp.toString(), JunctionTypeUp, fields);
            else if(truthsetValue.FieldName.equals(RnaFusionComparer.Fields.JuncTypeDown.toString()))
                addStringValue(RnaFusionComparer.Fields.JuncTypeDown.toString(), JunctionTypeDown, fields);
        }
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
    }

    @Override
    public String key()
    {
        return format("%s %s:%d-%s:%d", Name, ChromosomeUp, PositionUp, ChromosomeDown, PositionDown);
    }

    @Override
    public boolean reportable()
    {
        return KnownType != KnownFusionType.NONE;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaFusionData otherData = (RnaFusionData)other;

        if(!otherData.Name.equals(Name))
            return false;

        if(!otherData.ChromosomeUp.equals(ChromosomeUp))
            return false;

        if(!otherData.ChromosomeDown.equals(ChromosomeDown))
            return false;

        if(otherData.PositionUp != PositionUp)
            return false;

        return otherData.PositionDown == PositionDown;
    }
}
