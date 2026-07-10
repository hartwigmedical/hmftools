package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.common.CategoryType.V_CHORD;

import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record VChordData(VChordPrediction VChord) implements ComparableItem
{
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
