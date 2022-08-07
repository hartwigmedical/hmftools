package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.Category.CUPPA;
import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.Category.GERMLINE_DELETION;
import static com.hartwig.hmftools.compar.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.MismatchType.REF_ONLY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.cuppa.CuppaComparer;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.linx.DisruptionComparer;
import com.hartwig.hmftools.compar.linx.FusionComparer;
import com.hartwig.hmftools.compar.purple.GermlineDeletionComparer;
import com.hartwig.hmftools.compar.purple.PurityComparer;
import com.hartwig.hmftools.compar.somatic.GermlineVariantComparer;
import com.hartwig.hmftools.compar.somatic.SomaticVariantComparer;

public class CommonUtils
{
    public static final String DATA_DELIM = ",";
    public static final String ITEM_DELIM = ";";
    public static final String SUB_ITEM_DELIM = "=";

    public static final String FLD_REPORTED = "Reported";

    public static List<ItemComparer> buildComparers(final ComparConfig config)
    {
        List<ItemComparer> comparators = Lists.newArrayList();

        if(config.Categories.containsKey(PURITY))
            comparators.add(new PurityComparer(config));

        if(config.Categories.containsKey(DRIVER))
            comparators.add(new DriverComparer(config));

        if(config.Categories.containsKey(GERMLINE_DELETION))
            comparators.add(new GermlineDeletionComparer(config));

        if(config.Categories.containsKey(FUSION))
            comparators.add(new FusionComparer(config));

        if(config.Categories.containsKey(DISRUPTION))
            comparators.add(new DisruptionComparer(config));

        if(config.Categories.containsKey(SOMATIC_VARIANT))
            comparators.add(new SomaticVariantComparer(config));

        if(config.Categories.containsKey(GERMLINE_VARIANT))
            comparators.add(new GermlineVariantComparer(config));

        if(config.Categories.containsKey(CUPPA))
            comparators.add(new CuppaComparer(config));

        comparators.forEach(x -> x.registerThresholds(config.Thresholds));

        return comparators;
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
            for(int j = i + 1; j < config.SourceNames.size(); ++j)
            {
                CommonUtils.compareItems(mismatches, matchLevel, config.Thresholds, sourceItems.get(i), sourceItems.get(j));
            }
        }
    }

    public static void compareItems(
            final List<Mismatch> mismatches, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final List<ComparableItem> items1, final List<ComparableItem> items2)
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
                        Mismatch mismatch = item1.findMismatch(item2, matchLevel, thresholds);

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
                .forEach(x -> mismatches.add(new Mismatch(null, x, NEW_ONLY, emptyDiffs)));
    }
}
