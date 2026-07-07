package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.linx.BreakendData;

public class BreakendsField implements Field
{
    public final String name;
    public final Function<ComparableItem, List<BreakendData>> extractValue;
    public final boolean isCompared;

    public BreakendsField(final String name, final Function<ComparableItem, List<BreakendData>> extractValue, final boolean isCompared)
    {
        this.name = name;
        this.extractValue = extractValue;
        this.isCompared = isCompared;
    }

    @Override
    public String name()
    {
        return name;
    }

    @Override
    public boolean isCompared()
    {
        return isCompared;
    }

    @Override
    public String displayValue(final ComparableItem item)
    {
        return extractValue.apply(item).stream()
                .map(d -> d.fullStr(true))
                .collect(Collectors.joining(DISPLAY_VALUE_DELIMITER));
    }

    @Override
    public boolean hasDiff(final ComparableItem oldItem, final ComparableItem newItem)
    {
        return !determineDiffs(oldItem, newItem).isEmpty();
    }

    @Override
    public List<String> determineDiffs(ComparableItem oldItem, ComparableItem newItem)
    {
        final List<String> diffs = Lists.newArrayList();

        // compare each breakend and record differences
        List<BreakendData> oldBreakends = new ArrayList<>(extractValue.apply(oldItem));
        List<BreakendData> newBreakends = new ArrayList<>(extractValue.apply(newItem));

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
        return diffs;
    }

    private BreakendData findMatchingBreakend(final BreakendData breakend, final List<BreakendData> otherBreakends)
    {
        return otherBreakends.stream().filter(x -> x.matches(breakend)).findFirst().orElse(null);
    }
}
