package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.isofox.RnaGeneDataComparer.checkOldIsofoxFilename;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumMap;
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
        SplitFrags,
        Filter;
    }

    private final Map<SourceType,Map<String,List<RnaFusionData>>> mUnfilteredFusions;

    public RnaFusionComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mUnfilteredFusions = new EnumMap<>(SourceType.class);

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

        mFields.add(new FieldInfo(
                Fields.Filter.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Filter.toString()), null));
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
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

        loadUnfilteredFusions(sampleId, fileSources);

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(RnaFusion fusion : fusions)
        {
            RnaFusion passingFusion = ImmutableRnaFusion.builder().from(fusion).filter(PASS_FILTER).build();
            comparableItems.add(convertFusion(passingFusion, fileSources.Source));
        }

        return comparableItems;
    }

    private void loadUnfilteredFusions(final String sampleId, final PipelineSourcePaths fileSources)
    {
        mUnfilteredFusions.remove(fileSources.Source);

        String filename = RnaFusionFile.generateUnfilteredFilename(fileSources.Isofox, sampleId);
        filename = checkOldIsofoxFilename(filename);

        if(!Files.exists(Paths.get(filename)))
            return;

        List<RnaFusion> fusions;
        try
        {
            fusions = RnaFusionFile.read(filename);
        }
        catch(RuntimeException e)
        {
            CMP_LOGGER.warn("sample({}) failed to parse Isofox unfiltered fusion data ({})", sampleId, e.getClass().getSimpleName());
            return;
        }

        if(fusions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox unfiltered fusion data", sampleId);
            return;
        }

        List<RnaFusion> nonPassingFusions = fusions.stream().filter(x -> !x.filter().equals(PASS_FILTER)).collect(Collectors.toList());

        List<RnaFusionData> comparableFusions = nonPassingFusions.stream()
                .map(x -> convertFusion(x, fileSources.Source)).collect(Collectors.toList());

        mUnfilteredFusions.put(fileSources.Source, fusionsByName(comparableFusions));

        CMP_LOGGER.debug("sample({}) loaded {} {} non-passing fusions", sampleId, fileSources.Source, nonPassingFusions.size());
    }

    private RnaFusionData convertFusion(final RnaFusion fusion, final SourceType sourceType)
    {
        if(!mConfig.RequiresLiftover || sourceType != SourceType.OLD)
            return RnaFusionData.from(fusion, mFields);

        BasePosition liftedPositionUp = CommonUtils.determineComparisonGenomePosition(
                fusion.chromosomeUp(), fusion.positionUp(), sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);

        BasePosition liftedPositionDown = CommonUtils.determineComparisonGenomePosition(
                fusion.chromosomeDown(), fusion.positionDown(), sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);

        RnaFusion liftedFusion = ImmutableRnaFusion.builder().from(fusion)
                .chromosomeUp(liftedPositionUp.Chromosome)
                .chromosomeDown(liftedPositionDown.Chromosome)
                .positionUp(liftedPositionUp.Position)
                .positionDown(liftedPositionDown.Position)
                .build();

        return RnaFusionData.from(liftedFusion, mFields);
    }

    @Override
    public List<ComparableItem> loadFromTruthset(final Map<String,List<TruthsetValue>> valuesByKey)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(List<TruthsetValue> truthsetValues : valuesByKey.values())
        {
            comparableItems.add(RnaFusionData.fromTruthset(truthsetValues, mFields));
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
            checkUnfilteredFusions(mismatches);
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
                    oldFusionStr.add(fusionInfo(oldFusion, SourceType.NEW));
                }
                else
                {
                    Mismatch mismatch = oldFusion.findMismatch(this, newFusion, matchLevel, includeMatches, includesTruthset);

                    if(mismatch != null)
                        mismatches.add(mismatch);
                }
            }

            for(RnaFusionData newFusion : newFusions)
            {
                RnaFusionData oldFusion = oldFusions.stream().filter(x -> x.isMatched(newFusion)).findFirst().orElse(null);

                if(oldFusion == null)
                {
                    newFusionStr.add(fusionInfo(newFusion, SourceType.OLD));
                }
            }

            if(oldFusionStr.length() > 0 || newFusionStr.length() > 0)
            {
                ComparableItem oldItem = oldFusions.get(0);
                ComparableItem newItem = newFusions.get(0);
                String diffStr = format("fusions(%s/%s)", oldFusionStr, newFusionStr);
                Mismatch mismatch = new Mismatch(oldItem, newItem, VALUE, List.of(diffStr));
                mismatches.add(mismatch);
            }
        }

        // and add unmatched fusions
        addUnmatchedFusions(fusionsByNameOld, matchedByName, mismatches, matchLevel, false, oldTruthsetSourced);
        addUnmatchedFusions(fusionsByNameNew, matchedByName, mismatches, matchLevel, true, newTruthsetSourced);

        checkUnfilteredFusions(mismatches);
    }

    private void checkUnfilteredFusions(final List<Mismatch> mismatches)
    {
        for(int i = 0; i < mismatches.size(); ++i)
        {
            Mismatch mismatch = mismatches.get(i);

            boolean oldOnly = mismatch.Type == OLD_ONLY && mismatch.NewItem == null;
            boolean newOnly = mismatch.Type == NEW_ONLY && mismatch.OldItem == null;

            if(!oldOnly && !newOnly)
                continue;

            RnaFusionData fusion = (RnaFusionData)mismatch.nonNullItem();
            RnaFusionData unfilteredFusion = findUnfilteredFusion(fusion, oldOnly ? SourceType.NEW : SourceType.OLD);

            if(unfilteredFusion == null)
                continue;

            ComparableItem oldItem = oldOnly ? fusion : unfilteredFusion;
            ComparableItem newItem = oldOnly ? unfilteredFusion : fusion;
            List<String> diffs = CommonUtils.findDiffs(this, oldItem, newItem);
            mismatches.set(i, new Mismatch(oldItem, newItem, mismatch.Type, diffs));
        }
    }

    private RnaFusionData findUnfilteredFusion(final RnaFusionData fusion, final SourceType sourceType)
    {
        Map<String,List<RnaFusionData>> unfilteredFusions = mUnfilteredFusions.get(sourceType);

        if(unfilteredFusions == null)
            return null;

        return unfilteredFusions.getOrDefault(fusion.Name, Collections.emptyList()).stream()
                .filter(x -> x.isMatched(fusion)).findFirst().orElse(null);
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
                if(matchLevel == REPORTABLE && !fusion.reportable() && !truthsetSourced)
                    continue;

                Mismatch mismatch = isNew ?
                        new Mismatch(null, fusion, NEW_ONLY, Collections.emptyList()) :
                        new Mismatch(fusion, null, OLD_ONLY, Collections.emptyList());

                mismatches.add(mismatch);
            }
        }
    }

    private String fusionInfo(final RnaFusionData fusionData, final SourceType otherSource)
    {
        RnaFusionData unfilteredFusion = findUnfilteredFusion(fusionData, otherSource);

        if(unfilteredFusion == null)
            return fusionInfo(fusionData);

        return format("%s unfiltered(%s filter=%s)", fusionInfo(fusionData), fusionInfo(unfilteredFusion), unfilteredFusion.Filter);
    }

    private static String fusionInfo(final RnaFusionData fusionData)
    {
        return format("coords(%d-%d frags=%d)", fusionData.PositionUp, fusionData.PositionDown, fusionData.SplitFragments);
    }

    private static Map<String,List<RnaFusionData>> fusionsByName(final List<? extends ComparableItem> comparableItems)
    {
        Map<String,List<RnaFusionData>> fusionsByName = Maps.newHashMap();

        for(ComparableItem item : comparableItems)
        {
            RnaFusionData fusionData = (RnaFusionData)item;

            fusionsByName.computeIfAbsent(fusionData.Name, x -> Lists.newArrayList()).add(fusionData);
        }

        for(List<RnaFusionData> fusions : fusionsByName.values())
        {
            fusions.sort(Comparator.comparingInt((RnaFusionData x) -> x.SplitFragments).reversed());
        }

        return fusionsByName;
    }

}
