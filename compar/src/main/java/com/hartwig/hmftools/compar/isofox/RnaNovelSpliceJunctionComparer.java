package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
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
        String filename = determineFileName(sampleId, fileSources);
        List<NovelSpliceJunction> junctions = NovelSpliceJunctionFile.read(filename);
        if(junctions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Gene data", sampleId);
            return null;
        }

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(NovelSpliceJunction junction : junctions)
        {
            BasePosition comparisonPositionStart = CommonUtils.determineComparisonGenomePosition(
                    junction.chromosome(), junction.junctionStart(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            BasePosition comparisonPositionEnd = CommonUtils.determineComparisonGenomePosition(
                    junction.chromosome(), junction.junctionEnd(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            comparableItems.add(new RnaNovelSpliceJunctionData(junction, comparisonPositionStart, comparisonPositionEnd, mFields));
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final PipelineSourcePaths fileSources)
    {
        String filename = NovelSpliceJunctionFile.generateFilename(fileSources.Isofox, sampleId);
        String oldFilename = filename.replace(".tsv", ".csv");

        if(!Files.exists(Paths.get(filename)) && Files.exists(Paths.get(oldFilename)))
        {
            return oldFilename;
        }
        else
        {
            return filename;
        }
    }
}
