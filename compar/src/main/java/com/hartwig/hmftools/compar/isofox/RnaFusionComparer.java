package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.CategoryType.RNA_FUSION;
import static com.hartwig.hmftools.compar.isofox.RnaGeneDataComparer.checkOldIsofoxFilename;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
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
import com.hartwig.hmftools.compar.purple.PurityData;

public class RnaFusionComparer extends ItemComparer
{
    protected enum Fields
    {
        KnownFusionType,
        JuncTypeUp,
        JuncTypeDown,
        SplitFrags;
    }

    public RnaFusionComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.KnownFusionType.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.KnownFusionType.toString()), null));

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
        return false;
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

                comparableItems.add(new RnaFusionData(adjustedFusion, mFields));
            }
            else
            {
                comparableItems.add(new RnaFusionData(fusion, mFields));
           }
        }

        return comparableItems;
    }

    public List<ComparableItem> loadFromTruthset(final Map<String,List<TruthsetValue>> valuesByKey)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        // TODO
        
        return comparableItems;
    }
}
