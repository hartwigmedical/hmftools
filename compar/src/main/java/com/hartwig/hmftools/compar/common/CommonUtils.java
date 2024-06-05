package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.compar.common.Category.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_REF;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.REF_ONLY;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.chord.ChordComparer;
import com.hartwig.hmftools.compar.cuppa.CuppaComparer;
import com.hartwig.hmftools.compar.driver.DriverComparer;
import com.hartwig.hmftools.compar.lilac.LilacComparer;
import com.hartwig.hmftools.compar.linx.DisruptionComparer;
import com.hartwig.hmftools.compar.linx.FusionComparer;
import com.hartwig.hmftools.compar.linx.GermlineSvComparer;
import com.hartwig.hmftools.compar.purple.CopyNumberComparer;
import com.hartwig.hmftools.compar.purple.GeneCopyNumberComparer;
import com.hartwig.hmftools.compar.purple.GermlineDeletionComparer;
import com.hartwig.hmftools.compar.purple.PurityComparer;
import com.hartwig.hmftools.compar.mutation.GermlineVariantComparer;
import com.hartwig.hmftools.compar.mutation.SomaticVariantComparer;

public class CommonUtils
{
    public static final String SUB_ITEM_DELIM = "=";

    public static final String FLD_REPORTED = "Reported";
    public static final String FLD_QUAL = "Qual";

    public static List<ItemComparer> buildComparers(final ComparConfig config)
    {
        List<ItemComparer> comparers = Lists.newArrayList();

        // load in a predictable order irrespective of config
        for(Category category : Category.values())
        {
            if(!config.Categories.containsKey(category))
                continue;

            MatchLevel matchLevel = config.Categories.get(category);

            ItemComparer comparer = createComparer(category, config);

            if(matchLevel == REPORTABLE && !comparer.hasReportable())
                continue;

            comparer.registerThresholds(config.Thresholds);
            comparers.add(comparer);
        }

        if(config.runCopyNumberGeneComparer())
        {
            // use the gene copy number comparer to assist with CN-event drivers
            ItemComparer comparer = createComparer(GENE_COPY_NUMBER, config);
            comparer.registerThresholds(config.Thresholds);
            comparers.add(comparer);
        }

        return comparers;
    }

    private static ItemComparer createComparer(final Category category, final ComparConfig config)
    {
        switch(category)
        {
            case PURITY:
                return new PurityComparer(config);

            case DRIVER:
                return new DriverComparer(config);

            case COPY_NUMBER:
                return new CopyNumberComparer(config);

            case GENE_COPY_NUMBER:
                return new GeneCopyNumberComparer(config);

            case GERMLINE_DELETION:
                return new GermlineDeletionComparer(config);

            case FUSION:
                return new FusionComparer(config);

            case DISRUPTION:
                return new DisruptionComparer(config);

            case SOMATIC_VARIANT:
                return new SomaticVariantComparer(config);

            case GERMLINE_VARIANT:
                return new GermlineVariantComparer(config);

            case CUPPA:
                return new CuppaComparer(config);

            case CHORD:
                return new ChordComparer(config);

            case LILAC:
                return new LilacComparer(config);

            case GERMLINE_SV:
                return new GermlineSvComparer(config);

            default:
                return null;
        }
    }

    public static boolean processSample(
            final ItemComparer comparer, final ComparConfig config, final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = config.Categories.get(comparer.category());

        Map<String,List<ComparableItem>> sourceItems = Maps.newHashMap();

        for(String sourceName : config.SourceNames)
        {
            String sourceSampleId = config.sourceSampleId(sourceName, sampleId);
            List<ComparableItem> items = null;

            if(!config.DbConnections.isEmpty())
            {
                items = comparer.loadFromDb(sourceSampleId, config.DbConnections.get(sourceName), sourceName);
            }
            else
            {
                FileSources fileSources = FileSources.sampleInstance(config.FileSources.get(sourceName), sourceSampleId);
                items = comparer.loadFromFile(sourceSampleId, fileSources);
            }

            if(items != null)
            {
                sourceItems.put(sourceName, items);
            }
        }

        if(sourceItems.containsKey(REF_SOURCE) && sourceItems.containsKey(NEW_SOURCE))
        {
            // previously support comparisons for N sources but now can only be 2 as controlled by config
            CommonUtils.compareItems(mismatches, matchLevel, config.Thresholds, sourceItems.get(REF_SOURCE), sourceItems.get(NEW_SOURCE));
            return true;
        }

        InvalidDataItem invalidDataItem = new InvalidDataItem(comparer.category());

        if(!sourceItems.containsKey(REF_SOURCE) && !sourceItems.containsKey(NEW_SOURCE))
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, Collections.EMPTY_LIST));
        else if(!sourceItems.containsKey(REF_SOURCE))
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_REF, Collections.EMPTY_LIST));
        else if(!sourceItems.containsKey(NEW_SOURCE))
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, Collections.EMPTY_LIST));

        return false;
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

    public static BasePosition determineComparisonGenomePosition(
            final String chromosome, final int position, final String fileSource,
            final boolean requiresLiftover, final GenomeLiftoverCache liftoverCache)
    {
        if(requiresLiftover && liftoverCache.hasMappings())
        {
            RefGenomeVersion destVersion = fileSource.equals(REF_SOURCE) ? V38 : V37;
            int newPosition = liftoverCache.convertPosition(chromosome, position, destVersion);

            if(newPosition != UNMAPPED_POSITION)
            {
                return new BasePosition(destVersion.versionedChromosome(chromosome), newPosition);
            }
        }

        return new BasePosition(chromosome, position);
    }
}
