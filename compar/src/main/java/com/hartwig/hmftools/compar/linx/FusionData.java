package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;

import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class FusionData implements ComparableItem
{
    public final LinxFusion Fusion;
    public final String GeneMappedName;
    public final BreakendData BreakendFive;
    public final BreakendData BreakendThree;

    public FusionData(final LinxFusion fusion, final String geneMappedName, final BreakendData breakendFive, final BreakendData breakendThree)
    {
        Fusion = fusion;
        GeneMappedName = geneMappedName;
        BreakendFive = breakendFive;
        BreakendThree = breakendThree;
    }

    @Override
    public CategoryType category() { return FUSION; }

    @Override
    public String key()
    {
        return format("%s_%s", Fusion.name(), Fusion.reportedType());
    }

    @Override
    public boolean reportable() {
        return Fusion.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        return otherFusion.GeneMappedName.equals(GeneMappedName);
    }
}
