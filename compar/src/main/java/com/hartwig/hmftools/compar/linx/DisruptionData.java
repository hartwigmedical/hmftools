package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.field.FieldInfo.findField;

import java.util.List;

import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class DisruptionData extends ComparableItem
{
    public final String GeneName;
    public final List<BreakendData> Breakends;
    private final CategoryType mSubCategory;

    public DisruptionData(
            final CategoryType category, final String geneName, final List<BreakendData> breakends, final List<FieldInfo> fields)
    {
        mSubCategory = category;
        GeneName = geneName;
        Breakends = breakends;

        addBoolValue(DisruptionComparer.Fields.Reported.toString(), reportable(), fields);

        mValues.put(
                DisruptionComparer.Fields.BreakendInfo.toString(),
                new BreakendsFieldValue(findField(DisruptionComparer.Fields.BreakendInfo.toString(), fields), breakends));
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
