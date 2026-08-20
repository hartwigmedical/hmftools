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

    public RnaFusionData(final RnaFusion rnaFusion, final List<FieldInfo> fields)
    {
        Fusion = rnaFusion;

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

        if(!otherData.Fusion.chromosomeUp().equals(Fusion.chromosomeUp()))
            return false;

        if(!otherData.Fusion.chromosomeDown().equals(Fusion.chromosomeDown()))
            return false;

        if(otherData.Fusion.positionUp() != Fusion.positionUp())
            return false;

        return otherData.Fusion.positionDown() == Fusion.positionDown();
    }
}
