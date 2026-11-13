package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.compar.linx.BreakendData;

public class BreakendsFieldValue extends FieldValue
{
    private final List<BreakendData> mBreakends;

    public BreakendsFieldValue(final FieldInfo field, final List<BreakendData> breakends)
    {
        super(field);
        mBreakends = breakends;
    }

    @Override
    public boolean hasDifference(final FieldValue other)
    {
        BreakendsFieldValue otherValue = (BreakendsFieldValue)other;

        List<String> diffs = Lists.newArrayList();
        determineDiffs(mBreakends, otherValue.mBreakends, diffs);
        return !diffs.isEmpty();
    }

    @Override
    public String toString() { return format("%s=%s", name(), displayValue()); }

    @Override
    public String displayValue()
    {
        return mBreakends.stream().map(x -> x.fullStr(true)).collect(Collectors.joining(ITEM_DELIM));
    }

    public void addDiffInfo(final FieldValue oldValue, final FieldValue newValue, final List<String> diffs)
    {
        BreakendsFieldValue otherValue = (BreakendsFieldValue)newValue;
        determineDiffs(mBreakends, otherValue.mBreakends, diffs);
    }

    public void determineDiffs(final List<BreakendData> olds, final List<BreakendData> news, final List<String> diffs)
    {
        List<BreakendData> oldBreakends = Lists.newArrayList(olds);
        List<BreakendData> newBreakends = Lists.newArrayList(news);

        int index = 0;
        while(index < oldBreakends.size())
        {
            BreakendData oldBreakendData = oldBreakends.get(index);

            BreakendData newBreakendData = findMatchingBreakend(oldBreakendData, newBreakends);

            if(newBreakendData != null)
            {
                LinxBreakend oldBreakend = oldBreakendData.Breakend;
                LinxBreakend newBreakend = newBreakendData.Breakend;

                if(oldBreakend.regionType() != newBreakend.regionType()
                        || oldBreakend.codingType() != newBreakend.codingType()
                        || oldBreakend.nextSpliceExonRank() != newBreakend.nextSpliceExonRank())
                {
                    diffs.add(format("%s(%s transcript %s/%s)",
                            name(), oldBreakendData.svInfoStr(), oldBreakendData.transcriptStr(), newBreakendData.transcriptStr()));
                }

                if(oldBreakend.reportedStatus() != newBreakend.reportedStatus())
                {
                    diffs.add(format("%s(%s reported %s/%s)",
                            name(), oldBreakendData.svInfoStr(), oldBreakend.reportedStatus(), newBreakend.reportedStatus()));
                }

                oldBreakends.remove(index);
                newBreakends.remove(newBreakendData);
            }
            else
            {
                // record an unmatched breakend or SV
                diffs.add(format("%s(unmatched SV %s/)", name(), oldBreakendData.svInfoStr()));

                ++index;
            }
        }

        for(BreakendData otherBreakendData : newBreakends)
        {
            // record an unmatched breakend or SV on the other side
            diffs.add(format("%s(unmatched SV /%s)", name(), otherBreakendData.svInfoStr()));
        }
    }

    private BreakendData findMatchingBreakend(final BreakendData breakend, final List<BreakendData> otherBreakends)
    {
        return otherBreakends.stream().filter(x -> x.matches(breakend)).findFirst().orElse(null);
    }
}
