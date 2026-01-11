package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class DisruptionData implements ComparableItem
{
    public final String GeneName;
    public final List<BreakendData> Breakends;
    private final CategoryType mSubCategory;

    protected static final String FLD_BREAKEND_INFO = "BreakendInfo";
    protected static final String FLD_UNMATCHED_SV = "unmatchedSv";

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
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();

        values.add(String.format("%s", reportable()));

        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        for(BreakendData breakendData : Breakends)
        {
            sj.add(breakendData.fullStr());
        }

        values.add(sj.toString());
        return values;
    }

    @Override
    public boolean reportable() { return Breakends.stream().anyMatch(x -> x.Breakend.reportedStatus() == ReportedStatus.REPORTED); }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherDisruptionData = (DisruptionData)other;

        return GeneName.equals(otherDisruptionData.GeneName);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final DisruptionData otherDisruptionData = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        // compare each breakend and record differences
        List<BreakendData> breakends = Lists.newArrayList(Breakends);
        List<BreakendData> otherBreakends = Lists.newArrayList(otherDisruptionData.Breakends);

        int index = 0;
        while(index < breakends.size())
        {
            BreakendData breakendData = breakends.get(index);

            BreakendData otherBreakendData = findMatchingBreakend(breakendData, otherBreakends);

            if(otherBreakendData != null)
            {
                LinxBreakend breakend = breakendData.Breakend;
                LinxBreakend otherBreakend = otherBreakendData.Breakend;

                if(breakend.regionType() != otherBreakend.regionType()
                || breakend.codingType() != otherBreakend.codingType()
                || breakend.nextSpliceExonRank() != otherBreakend.nextSpliceExonRank())
                {
                    diffs.add(format("breakend(%s transcript %s/%s)",
                            breakendData.svInfoStr(), breakendData.transcriptStr(), otherBreakendData.transcriptStr()));
                }

                if(breakend.reportedStatus() != otherBreakend.reportedStatus())
                {
                    diffs.add(format("breakend(%s reported %s/%s)",
                            breakendData.svInfoStr(), breakend.reportedStatus(), otherBreakend.reportedStatus()));
                }

                breakends.remove(index);
                otherBreakends.remove(otherBreakendData);
            }
            else
            {
                // record an unmatched breakend or SV
                diffs.add(format(FLD_UNMATCHED_SV + "(%s/)", breakendData.svInfoStr()));

                ++index;
            }
        }

        for(BreakendData otherBreakendData : otherBreakends)
        {
            // record an unmatched breakend or SV on the other side
            diffs.add(format(FLD_UNMATCHED_SV + "(/%s)", otherBreakendData.svInfoStr()));
        }

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    public BreakendData findMatchingBreakend(final BreakendData breakend, final List<BreakendData> otherBreakends)
    {
        return otherBreakends.stream().filter(x -> x.matches(breakend)).findFirst().orElse(null);
    }
}
