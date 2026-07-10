package com.hartwig.hmftools.compar.teal;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.TELOMERE_LENGTH;

import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class TealData implements ComparableItem
{
    public final TelomereLength TelomereLength;

    public TealData(final TelomereLength telomereLength)
    {
        TelomereLength = telomereLength;
    }

    @Override
    public CategoryType category() { return TELOMERE_LENGTH; }

    @Override
    public String key()
    {
        return String.format("type(%s)", TelomereLength.type());
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final TelomereLength otherTelomereLength = ((TealData) other).TelomereLength;

        return TelomereLength.type().equals(otherTelomereLength.type());
    }

    public String toString()
    {
        return format("type(%s) telomere length(%.3f)", TelomereLength.type(), TelomereLength.finalTelomereLength());
    }
}
