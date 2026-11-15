package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.common.CategoryType.V_CHORD;

import java.util.List;

import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class VChordData extends ComparableItem
{
    public final VChordPrediction Prediction;

    VChordData(final VChordPrediction prediction, final List<FieldInfo> fields)
    {
        Prediction = prediction;

        addDoubleValue(VChordComparer.Fields.BreastCancerHrdScore.toString(), prediction.breastCancerHrdScore(), fields);
        addDoubleValue(VChordComparer.Fields.OvarianCancerHrdScore.toString(), prediction.ovarianCancerHrdScore(), fields);
        addDoubleValue(VChordComparer.Fields.OtherCancerScore.toString(), prediction.otherCancerScore(), fields);
    }

    @Override
    public CategoryType category()
    {
        return V_CHORD;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
