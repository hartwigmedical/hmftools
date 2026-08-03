package com.hartwig.hmftools.compar.isofox;

import java.util.List;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaFusionData extends ComparableItem
{
    public final RnaFusion Fusion;
    public final BasePosition ComparisonPositionUp;
    public final BasePosition ComparisonPositionDown;

    public RnaFusionData(
            final RnaFusion rnaFusion, final BasePosition comparisonPositionUp, final BasePosition comparisonPositionDown,
            final List<FieldInfo> fields)
    {
        Fusion = rnaFusion;
        ComparisonPositionUp = comparisonPositionUp;
        ComparisonPositionDown = comparisonPositionDown;

        addStringValue(RnaFusionComparer.Fields.KnownFusionType.toString(), rnaFusion.knownType().toString(), fields);
        addStringValue(RnaFusionComparer.Fields.JuncTypeUp.toString(), rnaFusion.junctionTypeUp(), fields);
        addStringValue(RnaFusionComparer.Fields.JuncTypeDown.toString(), rnaFusion.junctionTypeDown(), fields);
        addIntValue(RnaFusionComparer.Fields.SplitFrags.toString(), rnaFusion.splitFragments(), fields);
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%s:%d", Fusion.name(), Fusion.chromosomeUp(),
                Fusion.positionUp(), Fusion.chromosomeDown(), Fusion.positionDown());

        boolean upLifted = ComparisonPositionUp.Position != Fusion.positionUp()
                || !ComparisonPositionUp.Chromosome.equals(Fusion.chromosomeUp());

        boolean downLifted = ComparisonPositionDown.Position != Fusion.positionDown()
                || !ComparisonPositionDown.Chromosome.equals(Fusion.chromosomeDown());

        if(upLifted || downLifted)
            key += String.format(" liftover(%s-%s)", ComparisonPositionUp, ComparisonPositionDown);

        return key;
    }

    @Override
    public boolean reportable()
    {
        return Fusion.knownType() != KnownFusionType.NONE;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaFusionData otherData = (RnaFusionData)other;

        if(!otherData.Fusion.name().equals(Fusion.name())){
            return false;
        }
        if(!otherData.Fusion.chromosomeUp().equals(ComparisonPositionUp.Chromosome))
        {
            return false;
        }
        if(!otherData.Fusion.chromosomeDown().equals(ComparisonPositionDown.Chromosome))
        {
            return false;
        }
        if(otherData.Fusion.positionUp() != ComparisonPositionUp.Position)
        {
            return false;
        }
        return otherData.Fusion.positionDown() == ComparisonPositionDown.Position;
    }
}
