package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.FULL_MATCH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_ERROR;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_OLD;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.io.File;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.FieldValue;

public class CommonUtils
{
    public static boolean processSample(
            final ItemComparer comparer, final ComparConfig config, final String sampleId, final List<Mismatch> mismatches)
    {
        MatchLevel matchLevel = config.MatchingLevel;

        Map<SourceType,List<ComparableItem>> sourceItems = Maps.newHashMap();
        boolean includesTruthset = false;

        for(SourceData source : config.Sources)
        {
            String sourceSampleId = config.sourceSampleId(source.Type, sampleId);
            String sourceReferenceId = config.sourceReferenceId(source.Type, sampleId);
            List<ComparableItem> items = null;

            if(source.PipelinePaths != null)
            {
                PipelineSourcePaths fileSources = PipelineSourcePaths.sampleInstance(source.PipelinePaths, sourceSampleId, sourceReferenceId);
                items = comparer.loadFromFile(sourceSampleId, sourceReferenceId, fileSources);
            }
            else
            {
                includesTruthset = true;

                List<TruthsetValue> truthsetValues = source.Truthset.sampleTruthsetEntries(sampleId, comparer.category());

                if(truthsetValues == null || truthsetValues.isEmpty())
                    continue;

                // validate truthset fields against the comparer
                List<FieldInfo> fields = comparer.fieldsList();

                for(TruthsetValue truthsetValue : truthsetValues)
                {
                    if(fields.stream().noneMatch(x -> x.Name.equals(truthsetValue.FieldName)))
                    {
                        CMP_LOGGER.error("category({}) invalid truthset entry({})", comparer.category(), truthsetValue);
                        return false;
                    }
                }

                // group by items keys
                Map<String, List<TruthsetValue>> truthsetValuesByKey = truthsetValues.stream().collect(Collectors.groupingBy(x -> x.Key));

                items = comparer.loadFromTruthset(truthsetValuesByKey);
            }

            if(items != null)
            {
                sourceItems.put(source.Type, items);
            }
        }

        if(sourceItems.containsKey(SourceType.OLD) && sourceItems.containsKey(SourceType.NEW))
        {
            comparer.compareItems(
                    mismatches, matchLevel, config.IncludeMatches, includesTruthset,
                    sourceItems.get(SourceType.OLD), sourceItems.get(SourceType.NEW));

            return true;
        }

        InvalidDataItem invalidDataItem = new InvalidDataItem(comparer.category());

        if(!sourceItems.containsKey(SourceType.OLD) && !sourceItems.containsKey(SourceType.NEW))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_BOTH, Collections.EMPTY_LIST));
        }
        else if(!sourceItems.containsKey(SourceType.OLD))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_OLD, Collections.EMPTY_LIST));
        }
        else if(!sourceItems.containsKey(SourceType.NEW))
        {
            mismatches.add(new Mismatch(invalidDataItem, null, INVALID_NEW, Collections.EMPTY_LIST));
        }

        return false;
    }

    public static void compareItems(
            final ItemComparer comparer, final List<Mismatch> mismatches, final MatchLevel matchLevel, final boolean includeMatches,
            final boolean includesTruthset, final List<ComparableItem> oldItems, final List<ComparableItem> newItems)
    {
        boolean oldTruthsetSourced = comparer.config().isTruthsetSourced(SourceType.OLD);
        boolean newTruthsetSourced = comparer.config().isTruthsetSourced(SourceType.NEW);

        int index1 = 0;
        while(index1 < oldItems.size())
        {
            final ComparableItem item1 = oldItems.get(index1);

            boolean matched = false;

            int index2 = 0;
            while(index2 < newItems.size())
            {
                final ComparableItem item2 = newItems.get(index2);

                if(item1.matches(item2))
                {
                    oldItems.remove(index1);
                    newItems.remove(index2);
                    matched = true;

                    // skip checking for diffs if the items are not reportable
                    boolean checkMismatch;

                    if(matchLevel != REPORTABLE)
                    {
                        checkMismatch = true;
                    }
                    else if(item1.reportable() || item2.reportable())
                    {
                        checkMismatch = true;
                    }
                    else
                    {
                        checkMismatch = oldTruthsetSourced || newTruthsetSourced;
                    }

                    if(checkMismatch)
                    {
                        Mismatch mismatch = item1.findMismatch(comparer, item2, matchLevel, includeMatches, includesTruthset);

                        if(mismatch != null)
                        {
                            mismatches.add(mismatch);
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
            {
                ++index1;
            }
        }

        if(oldItems.isEmpty() && newItems.isEmpty())
        {
            return;
        }

        List<String> emptyDiffs = Lists.newArrayList();

        oldItems.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(x, null, OLD_ONLY, emptyDiffs)));

        newItems.stream().filter(x -> matchLevel != REPORTABLE || x.reportable())
                .forEach(x -> mismatches.add(new Mismatch(null, x, NEW_ONLY, emptyDiffs)));
    }

    public static List<String> findDiffs(final ItemComparer comparer, final ComparableItem oldItem, final ComparableItem newItem)
    {
        List<String> diffs = Lists.newArrayList();

        // find and compare fields present in both items
        List<FieldInfo> fields = comparer.fieldsList();
        Map<String,FieldValue> oldFieldValues = oldItem.fieldValues();
        Map<String,FieldValue> newFieldValues = newItem.fieldValues();

        for(FieldInfo field : fields)
        {
            if(!field.FieldCheck.IsCompared)
                continue;

            FieldValue oldValue = oldFieldValues.get(field.Name);
            FieldValue newValue = newFieldValues.get(field.Name);

            if(oldValue == null || newValue == null)
                continue;

            if(oldValue.hasDifference(newValue))
            {
                oldValue.addDiffInfo(oldValue, newValue, diffs);
            }
        }

        return diffs;
    }

    public static BasePosition determineComparisonGenomePosition(
            final String chromosome, final int position, final SourceType sourceType,
            final boolean requiresLiftover, final GenomeLiftoverCache liftoverCache)
    {
        if(requiresLiftover && liftoverCache.hasMappings())
        {
            RefGenomeVersion destVersion = sourceType == SourceType.OLD ? V38 : V37;
            int newPosition = liftoverCache.convertPosition(chromosome, position, destVersion);

            if(newPosition != UNMAPPED_POSITION)
            {
                return new BasePosition(destVersion.versionedChromosome(chromosome), newPosition);
            }
        }

        return new BasePosition(chromosome, position);
    }

    public static String determineComparisonChromosome(final String chromosome, final boolean requiresLiftover)
    {
        if(requiresLiftover)
        {
            return HumanChromosome.fromString(chromosome).name().substring(1);
        }
        else
        {
            return chromosome;
        }
    }

    public static boolean fileExists(final String filename)
    {
        return Files.exists(new File(filename).toPath());
    }

    public static Mismatch createMismatchFromDiffs(
            final ComparableItem oldItem, final ComparableItem newItem, final List<String> diffs,
            final MatchLevel matchLevel, final boolean includeMatches, final boolean includesTruthset)
    {
        boolean oldCountsAsCalled = countsAsCalled(oldItem, matchLevel) || includesTruthset;
        boolean newCountsAsCalled = countsAsCalled(newItem, matchLevel) || includesTruthset;

        if(!oldCountsAsCalled && !newCountsAsCalled)
        {
            // ignore unimportant differences
            return null;
        }
        else if(oldCountsAsCalled && newCountsAsCalled && diffs.isEmpty() && !includeMatches)
        {
            // ignore perfect matches when not including matches
            return null;
        }
        else
        {
            MismatchType mismatchType;
            if(oldCountsAsCalled && !newCountsAsCalled)
            {
                mismatchType = OLD_ONLY;
            }
            else if(!oldCountsAsCalled && newCountsAsCalled)
            {
                mismatchType = NEW_ONLY;
            }
            else if(oldCountsAsCalled && newCountsAsCalled && !diffs.isEmpty())
            {
                mismatchType = VALUE;
            }
            else if(oldCountsAsCalled && newCountsAsCalled && includeMatches)
            {
                mismatchType = FULL_MATCH;
            }
            else
            {
                // should be impossible due to earlier filters
                mismatchType = INVALID_ERROR;
            }

            return new Mismatch(oldItem, newItem, mismatchType, diffs);
        }
    }

    public static boolean countsAsCalled(final ComparableItem item, final MatchLevel matchLevel)
    {
        if(matchLevel == REPORTABLE)
        {
            return item.reportable();
        }
        else
        {
            return item.isPass();
        }
    }
}
