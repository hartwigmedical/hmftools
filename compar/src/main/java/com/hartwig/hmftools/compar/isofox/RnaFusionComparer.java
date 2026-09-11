package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.FULL_MATCH;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.isofox.RnaGeneDataComparer.checkOldIsofoxFilename;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaFusionComparer extends ItemComparer
{
    protected enum Fields
    {
        KnownType,
        JuncTypeUp,
        JuncTypeDown,
        SplitFrags;
    }

    public RnaFusionComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.KnownType.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.KnownType.toString()), null));

        mFields.add(new FieldInfo(
                Fields.JuncTypeUp.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.JuncTypeUp.toString()), null));

        mFields.add(new FieldInfo(
                Fields.JuncTypeDown.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.JuncTypeDown.toString()), null));

        mFields.add(new FieldInfo(
                Fields.SplitFrags.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SplitFrags.toString(), 5.0, 0.05),
                "%.2f"));
    }
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
    }

    @Override
    public boolean hasReportable()
    {
        return true;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        String filename = RnaFusionFile.generateFilename(fileSources.Isofox, sampleId);
        filename = checkOldIsofoxFilename(filename);

        List<RnaFusion> fusions = RnaFusionFile.read(filename);
        if(fusions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox fusion data", sampleId);
            return null;
        }

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(RnaFusion fusion : fusions)
        {
            if(mConfig.RequiresLiftover && fileSources.Source == SourceType.OLD)
            {
                BasePosition liftedPositionUp = CommonUtils.determineComparisonGenomePosition(
                        fusion.chromosomeUp(), fusion.positionUp(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                BasePosition liftedPositionDown = CommonUtils.determineComparisonGenomePosition(
                        fusion.chromosomeDown(), fusion.positionDown(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                RnaFusion adjustedFusion = ImmutableRnaFusion.builder().from(fusion)
                        .chromosomeUp(liftedPositionUp.Chromosome)
                        .chromosomeDown(liftedPositionDown.Chromosome)
                        .positionUp(liftedPositionUp.Position)
                        .positionDown(liftedPositionDown.Position)
                        .build();

                comparableItems.add(RnaFusionData.from(adjustedFusion, mFields));
            }
            else
            {
                comparableItems.add(RnaFusionData.from(fusion, mFields));
           }
        }

        return comparableItems;
    }

    @Override
    public List<ComparableItem> loadFromTruthset(final Map<String,List<TruthsetValue>> valuesByKey)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(List<TruthsetValue> truthsetValues : valuesByKey.values())
        {
            RnaFusionData fusionData = RnaFusionData.fromTruthset(truthsetValues, mFields);
            comparableItems.add(fusionData);
        }

        return comparableItems;
    }

    @Override
    public void compareItems(
            final List<Mismatch> mismatches, final MatchLevel matchLevel, final boolean includeMatches,
            final boolean includesTruthset, final List<ComparableItem> oldItems, final List<ComparableItem> newItems)
    {
        boolean oldTruthsetSourced = mConfig.isTruthsetSourced(SourceType.OLD);
        boolean newTruthsetSourced = mConfig.isTruthsetSourced(SourceType.NEW);

        if(!oldTruthsetSourced && !newTruthsetSourced)
        {
            CommonUtils.compareItems(this, mismatches, matchLevel, includeMatches, includesTruthset, oldItems, newItems);
            return;
        }

        // first organise fusion by name
        Map<String,List<RnaFusionData>> fusionsByNameOld = fusionsByName(oldItems);
        Map<String,List<RnaFusionData>> fusionsByNameNew = fusionsByName(newItems);

        Set<String> matchedByName = Sets.newHashSet();

        for(Map.Entry<String,List<RnaFusionData>> oldEntry : fusionsByNameOld.entrySet())
        {
            String fusionName = oldEntry.getKey();
            List<RnaFusionData> newFusions = fusionsByNameNew.get(fusionName);
            List<RnaFusionData> oldFusions = oldEntry.getValue();

            if(newFusions == null)
                continue;

            matchedByName.add(fusionName);

            StringJoiner oldFusionStr = new StringJoiner(ITEM_DELIM);
            StringJoiner newFusionStr = new StringJoiner(ITEM_DELIM);

            // match on name, then on coordinates
            // for truthset entries, if there is at least one match on coordinates then report other differences as value differences
            for(RnaFusionData oldFusion : oldFusions)
            {
                RnaFusionData newFusion = newFusions.stream().filter(x -> x.isMatched(oldFusion)).findFirst().orElse(null);

                if(newFusion == null)
                {
                    oldFusionStr.add(fusionInfo(oldFusion));
                }
                else
                {
                    mismatches.add(oldFusion.findMismatch(this, newFusion, matchLevel, includeMatches, includesTruthset));
                }
            }

            for(RnaFusionData newFusion : newFusions)
            {
                RnaFusionData oldFusion = oldFusions.stream().filter(x -> x.isMatched(newFusion)).findFirst().orElse(null);

                if(oldFusion == null)
                {
                    newFusionStr.add(fusionInfo(newFusion));
                }
            }

            if(oldFusionStr.length() > 0 || newFusionStr.length() > 0)
            {
                ComparableItem oldItem = !oldFusions.isEmpty() ? oldFusions.get(0) : null;
                ComparableItem newItem = !newFusions.isEmpty() ? newFusions.get(0) : null;
                String diffStr = format("fusions(%s/%s)", oldFusionStr, newFusionStr);
                Mismatch mismatch = new Mismatch(oldItem, newItem, VALUE, List.of(diffStr));
                mismatches.add(mismatch);
            }
        }

        // and add unmatched fusions
        addUnmatchedFusions(fusionsByNameOld, matchedByName, mismatches, matchLevel, false, oldTruthsetSourced);
        addUnmatchedFusions(fusionsByNameNew, matchedByName, mismatches, matchLevel, true, newTruthsetSourced);
    }

    private static void addUnmatchedFusions(
            final Map<String,List<RnaFusionData>> fusionsByName, final Set<String> matchedByName,
            final List<Mismatch> mismatches, final MatchLevel matchLevel, boolean isNew, boolean truthsetSourced)
    {
        for(Map.Entry<String,List<RnaFusionData>> entry : fusionsByName.entrySet())
        {
            String fusionName = entry.getKey();

            if(matchedByName.contains(fusionName))
                continue;

            for(RnaFusionData fusion : entry.getValue())
            {
                if((matchLevel == REPORTABLE && !fusion.reportable() && !truthsetSourced))
                    continue;

                Mismatch mismatch = isNew ?
                        new Mismatch(null, fusion, NEW_ONLY, Collections.emptyList()) :
                        new Mismatch(fusion, null, OLD_ONLY, Collections.emptyList());

                mismatches.add(mismatch);
            }
        }
    }

    private static String fusionInfo(final RnaFusionData fusionData)
    {
        return format("coords(%d-%d frags=%d)", fusionData.PositionUp, fusionData.PositionDown, fusionData.SplitFragments);
    }

    private static Map<String,List<RnaFusionData>> fusionsByName(final List<ComparableItem> comparableItems)
    {
        Map<String,List<RnaFusionData>> fusionsByName = Maps.newHashMap();

        for(ComparableItem item : comparableItems)
        {
            RnaFusionData fusionData = (RnaFusionData)item;

            List<RnaFusionData> fusions = fusionsByName.get(fusionData.Name);

            if(fusions == null)
            {
                fusions = Lists.newArrayList();
                fusionsByName.put(fusionData.Name, fusions);
            }

            fusions.add(fusionData);
        }

        for(List<RnaFusionData> fusions : fusionsByName.values())
        {
            Collections.sort(fusions, Comparator.comparingInt(x -> -x.SplitFragments));
        }

        return fusionsByName;
    }

}
