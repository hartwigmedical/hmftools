package com.hartwig.hmftools.compar.linx;

import java.util.List;

import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class DisruptionData implements ComparableItem
{
    public final String GeneName;
    public final List<BreakendData> Breakends;
    private final CategoryType mSubCategory;

    protected static final String FLD_BREAKEND = "Breakend";

    public DisruptionData(final CategoryType category, final String geneName, final List<BreakendData> breakends)
    {
        mSubCategory = category;
        GeneName = geneName;
        Breakends = breakends;
    }

    @Override
    public CategoryType category() { return mSubCategory; }

    @Override
    public String key()
    {
        return String.format("%s breakends(%d)", GeneName, Breakends.size());
    }

    @Override
    public boolean reportable() { return Breakends.stream().anyMatch(x -> x.Breakend.reportedStatus() == ReportedStatus.REPORTED); }

    @Override
    public String geneName() { return GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherDisruptionData = (DisruptionData)other;

        return GeneName.equals(otherDisruptionData.GeneName);
    }
}
