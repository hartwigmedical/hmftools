package com.hartwig.hmftools.compar;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.compar.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.Category.SV;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.PRESENCE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.linx.DisruptionComparer;
import com.hartwig.hmftools.compar.linx.FusionComparer;
import com.hartwig.hmftools.compar.linx.LinxSvComparer;
import com.hartwig.hmftools.compar.purple.CopyNumberComparer;
import com.hartwig.hmftools.compar.purple.PurityComparer;
import com.hartwig.hmftools.compar.somatic.SomaticVariantComparer;

public class CommonUtils
{
    public static final String DATA_DELIM = ",";
    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = "=";

    public static final double DEFAULT_DIFF = 0.1;
    public static final double DEFAULT_DIFF_PERC = 0.1;

    public static List<ItemComparer> buildComparers(final ComparConfig config)
    {
        List<ItemComparer> comparators = Lists.newArrayList();

        if(config.Categories.containsKey(PURITY))
            comparators.add(new PurityComparer(config));

        if(config.Categories.containsKey(COPY_NUMBER))
            comparators.add(new CopyNumberComparer(config));

        if(config.Categories.containsKey(DRIVER))
            comparators.add(new DriverComparer(config));

        if(config.Categories.containsKey(LINX_DATA) || config.Categories.containsKey(SV))
            comparators.add(new LinxSvComparer(config));

        if(config.Categories.containsKey(FUSION))
            comparators.add(new FusionComparer(config));

        if(config.Categories.containsKey(DISRUPTION))
            comparators.add(new DisruptionComparer(config));

        if(config.Categories.containsKey(SOMATIC_VARIANT))
            comparators.add(new SomaticVariantComparer(config));

        return comparators;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, double value1, double value2)
    {
        if(!diffValue(value1, value2))
            return false;

        diffs.add(String.format("%s(%.3f/%.3f)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, int value1, int value2)
    {
        if(value1 == value2)
            return false;

        diffs.add(String.format("%s(%d/%d)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, boolean value1, boolean value2)
    {
        if(value1 == value2)
            return false;

        diffs.add(String.format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, String value1, final String value2)
    {
        if(value1.equals(value2))
            return false;

        diffs.add(String.format("%s(%s/%s)", name, value1, value2));
        return true;
    }

    public static boolean diffValue(double value1, double value2)
    {
        return diffValue(value1, value2, DEFAULT_DIFF, DEFAULT_DIFF_PERC, true);
    }

    public static boolean diffValue(double value1, double value2, double maxAbsDiff, double maxRelDiff, boolean requireBothDiff)
    {
        // return TRUE if the values are DIFFERENT
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

    public static void processSample(
            final ItemComparer comparer, final ComparConfig config, final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = config.Categories.get(comparer.category());

        final List<List<ComparableItem>> sourceItems = Lists.newArrayList();

        for(String sourceName : config.SourceNames)
        {
            String sourceSampleId = config.sourceSampleId(sourceName, sampleId);

            if(!config.DbConnections.isEmpty())
                sourceItems.add(comparer.loadFromDb(sourceSampleId, config.DbConnections.get(sourceName)));
            else
                sourceItems.add(comparer.loadFromFile(sourceSampleId, config.FileSources.get(sourceName)));
        }

        for(int i = 0; i < config.SourceNames.size() - 1; ++i)
        {
            final String source1 = config.SourceNames.get(i);

            for(int j = i + 1; j < config.SourceNames.size(); ++j)
            {
                final String source2 = config.SourceNames.get(j);

                CommonUtils.compareItems(mismatches, matchLevel, source1, source2, sourceItems.get(i), sourceItems.get(j));
            }
        }
    }

    public static void compareItems(
            final List<Mismatch> mismatches, final MatchLevel matchLevel,
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
                    boolean eitherReportable = item2.reportable() || item2.reportable();

                    if(matchLevel != REPORTABLE || eitherReportable)
                    {
                        final List<String> diffs = item1.findDifferences(item2, matchLevel);

                        if(!diffs.isEmpty())
                        {
                            StringJoiner differencesStr = new StringJoiner(ITEM_DELIM);
                            diffs.forEach(x -> differencesStr.add(x));

                            mismatches.add(new Mismatch(
                                    item1.category(), VALUE, source1, source2, eitherReportable,
                                    item1.description(), differencesStr.toString()));
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
                .forEach(x -> mismatches.add(new Mismatch(x.category(), PRESENCE, source1, source2, x.reportable(), x.description(), "")));

        items2.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x.category(), PRESENCE, source2, source1, x.reportable(), x.description(), "")));
    }

}
