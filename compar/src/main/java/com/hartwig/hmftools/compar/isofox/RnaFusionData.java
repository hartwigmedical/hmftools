package com.hartwig.hmftools.compar.isofox;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record RnaFusionData(RnaFusion RnaFusion, BasePosition ComparisonPositionUp, BasePosition ComparisonPositionDown)
        implements ComparableItem
{
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
    }

    @Override
    public String key()
    {
        String key = String.format("%s %s:%d-%s:%d", RnaFusion.name(), RnaFusion.chromosomeUp(),
                RnaFusion.positionUp(), RnaFusion.chromosomeDown(), RnaFusion.positionDown());

        boolean upLifted = ComparisonPositionUp.Position != RnaFusion.positionUp()
                || !ComparisonPositionUp.Chromosome.equals(RnaFusion.chromosomeUp());

        boolean downLifted = ComparisonPositionDown.Position != RnaFusion.positionDown()
                || !ComparisonPositionDown.Chromosome.equals(RnaFusion.chromosomeDown());

        if(upLifted || downLifted)
            key += String.format(" liftover(%s-%s)", ComparisonPositionUp, ComparisonPositionDown);

        return key;
    }

    @Override
    public boolean reportable()
    {
        return RnaFusion.knownType() != KnownFusionType.NONE;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final RnaFusionData otherData = (RnaFusionData)other;

        if(!otherData.RnaFusion.name().equals(RnaFusion.name())){
            return false;
        }
        if(!otherData.RnaFusion.chromosomeUp().equals(ComparisonPositionUp.Chromosome))
        {
            return false;
        }
        if(!otherData.RnaFusion.chromosomeDown().equals(ComparisonPositionDown.Chromosome))
        {
            return false;
        }
        if(otherData.RnaFusion.positionUp() != ComparisonPositionUp.Position)
        {
            return false;
        }
        return otherData.RnaFusion.positionDown() == ComparisonPositionDown.Position;
    }
}
