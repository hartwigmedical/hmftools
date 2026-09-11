package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.isofox.RnaGeneDataComparer.checkOldIsofoxFilename;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaNovelSpliceJunctionComparer extends ItemComparer
{
    protected enum Fields
    {
        Type,
        FragmentCount,
        RegionStart,
        RegionEnd;
    }

    public RnaNovelSpliceJunctionComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Type.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Type.toString()), null));

        mFields.add(new FieldInfo(
                Fields.FragmentCount.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FragmentCount.toString(), 5.0, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.RegionStart.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.RegionStart.toString()), null));

        mFields.add(new FieldInfo(
                Fields.RegionEnd.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.RegionEnd.toString()), null));
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_NOVEL_SPLICE_JUNCTION;
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
        String filename = NovelSpliceJunctionFile.generateFilename(fileSources.Isofox, sampleId);
        filename = checkOldIsofoxFilename(filename);

        List<NovelSpliceJunction> junctions = NovelSpliceJunctionFile.read(filename);
        if(junctions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox alt splice junction data", sampleId);
            return null;
        }

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(NovelSpliceJunction junction : junctions)
        {
            if(mConfig.RequiresLiftover && fileSources.Source == SourceType.OLD)
            {
                BasePosition liftoverPositionStart = CommonUtils.determineComparisonGenomePosition(
                        junction.chromosome(), junction.junctionStart(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                BasePosition liftoverPositionEnd = CommonUtils.determineComparisonGenomePosition(
                        junction.chromosome(), junction.junctionEnd(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                NovelSpliceJunction adjustedJunction = ImmutableNovelSpliceJunction.builder().from(junction)
                        .chromosome(liftoverPositionStart.Chromosome)
                        .junctionStart(liftoverPositionStart.Position)
                        .junctionEnd(liftoverPositionEnd.Position)
                        .build();

                comparableItems.add(RnaNovelSpliceJunctionData.from(adjustedJunction, mFields));
            }
            else
            {
                comparableItems.add(RnaNovelSpliceJunctionData.from(junction, mFields));
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
            RnaNovelSpliceJunctionData junction = RnaNovelSpliceJunctionData.fromTruthset(truthsetValues, mFields);
            comparableItems.add(junction);
        }

        return comparableItems;
    }

}
