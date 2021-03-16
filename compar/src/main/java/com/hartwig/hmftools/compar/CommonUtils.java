package com.hartwig.hmftools.compar;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.PRESENCE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;
import java.util.StringJoiner;

public class CommonUtils
{
    public static final String DATA_DELIM = ",";
    public static final String ITEM_DELIM = ";";


    public static final double DEFAULT_DIFF = 0.1;
    public static final double DEFAULT_DIFF_PERC = 0.1;

    public static boolean diffValue(double value1, double value2)
    {
        return diffValue(value1, value2, DEFAULT_DIFF, DEFAULT_DIFF_PERC, true);
    }

    public static boolean diffValue(double value1, double value2, double maxAbsDiff, double maxRelDiff, boolean requireBothDiff)
    {
        if(value1 == 0 && value2 == 0)
            return false;

        double absDiff = abs(value1 - value2);
        boolean hasAbsDiff = absDiff > maxAbsDiff;
        boolean hasRelDiff = absDiff / max(value1, value2) > maxRelDiff;

        if(requireBothDiff)
            return hasAbsDiff && hasRelDiff;
        else
            return hasAbsDiff || hasRelDiff;
    }

    public static void compareItems(
            final String sampleId, final List<Mismatch> mismatches, final MatchLevel matchLevel,
            final String source1, final String source2, final List<ComparableItem> items1, final List<ComparableItem> items2)
    {
        int index1 = 0;
        while(index1 < items1.size())
        {
            final ComparableItem item1 = items1.get(index1);

            boolean matched = false;

            int index2 = 0;
            while(index2 < items2.size())
            {
                final ComparableItem item2 = items2.get(index2);

                if(item1.matches(item2))
                {
                    items1.remove(index1);
                    items2.remove(index2);
                    matched = true;

                    // skip checking for diffs if the items are not reportable
                    if(matchLevel != REPORTABLE || item2.reportable() || item2.reportable())
                    {
                        final List<String> diffs = item1.findDifferences(item2, matchLevel);

                        if(!diffs.isEmpty())
                        {
                            StringJoiner differencesStr = new StringJoiner(ITEM_DELIM);
                            diffs.forEach(x -> differencesStr.add(x));

                            mismatches.add(new Mismatch(
                                    item1.category(), VALUE, source1, source2, item1.description(), differencesStr.toString()));
                        }
                    }

                    break;
                }
                else
                {
                    ++index2;
                }
            }

            if(!matched)
                ++index1;
        }

        items1.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x.category(), PRESENCE, source1, source2, x.description(), "")));

        items2.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x.category(), PRESENCE, source2, source1, x.description(), "")));
    }

}
