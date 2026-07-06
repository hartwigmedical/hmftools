package com.hartwig.hmftools.compar.teal;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.TELOMERE_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class TealData implements ComparableItem
{
    public final TelomereLength TelomereLength;

    protected static final String FLD_TYPE = "Type";
    protected static final String FLD_TELOMERE_LENGTH = "TelomereLength";

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

    public static List<String> comparedFieldNames()
    {
        return List.of(FLD_TELOMERE_LENGTH);
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
