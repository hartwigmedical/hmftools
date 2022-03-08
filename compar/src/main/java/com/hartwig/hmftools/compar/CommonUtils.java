package com.hartwig.hmftools.compar;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.compar.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.Category.SV;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.MismatchType.REF_ONLY;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.linx.DisruptionComparer;
import com.hartwig.hmftools.compar.linx.FusionComparer;
import com.hartwig.hmftools.compar.linx.LinxSvComparer;
import com.hartwig.hmftools.compar.purple.CopyNumberComparer;
import com.hartwig.hmftools.compar.purple.PurityComparer;
import com.hartwig.hmftools.compar.somatic.GermlineVariantComparer;
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

        if(config.Categories.containsKey(GERMLINE_VARIANT))
            comparators.add(new GermlineVariantComparer(config));

        return comparators;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, double value1, double value2)
    {
        if(!hasSignificantDiff(value1, value2))
            return false;

        diffs.add(String.format("%s(%.3f/%.3f)", name, value1, value2));
        return true;
    }

    public static boolean checkDiff(final List<String> diffs, final String name, int value1, int value2)
    {
        if(!hasSignificantDiff(value1, value2))
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

    public static boolean hasSignificantDiff(double value1, double value2)
    {
        return hasSignificantDiff(value1, value2, DEFAULT_DIFF, DEFAULT_DIFF_PERC, true);
    }

    public static boolean hasSignificantDiff(double value1, double value2, double maxAbsDiff, double maxRelDiff, boolean requireBothDiff)
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
            {
                sourceItems.add(comparer.loadFromDb(sourceSampleId, config.DbConnections.get(sourceName)));
            }
            else
            {
                FileSources fileSources = config.FileSources.get(sourceName);
                sourceItems.add(comparer.loadFromFile(sourceSampleId, FileSources.sampleInstance(fileSources, sampleId)));
            }
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
                        // final List<String> diffs = item1.findDifferences(item2, matchLevel);
                        Mismatch mismatch = item1.findMismatch(item2, matchLevel);

                        if(mismatch != null)
                            mismatches.add(mismatch);
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

        if(items1.isEmpty() && items2.isEmpty())
            return;

        List<String> emptyDiffs = Lists.newArrayList();

        items1.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x, null, REF_ONLY, emptyDiffs)));

        items2.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x, null, NEW_ONLY, emptyDiffs)));
    }

    public static String filtersStr(final Set<String> filters)
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    public static String diffsStr(final List<String> diffs)
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        diffs.forEach(x -> sj.add(x));
        return sj.toString();
    }

}
