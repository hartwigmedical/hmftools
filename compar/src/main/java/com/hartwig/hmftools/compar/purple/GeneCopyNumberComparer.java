package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_PERC_THRESHOLD;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GeneCopyNumberComparer extends ItemComparer
{
    protected enum Fields
    {
        MinCopyNumber,
        MaxCopyNumber,
        DriverType,
        MinRegionStart,
        MinRegionEnd,
        TranscriptId,
        IsCanonical;
    }

    public GeneCopyNumberComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.MinCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MinCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MaxCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MaxCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.DriverType.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.DriverType.toString()), null));

        mFields.add(FieldInfo.displayOnly(Fields.MinRegionStart.toString(), null));
        mFields.add(FieldInfo.displayOnly(Fields.MinRegionEnd.toString(), null));
        mFields.add(FieldInfo.displayOnly(Fields.TranscriptId.toString(), null));
        mFields.add(FieldInfo.displayOnly(Fields.IsCanonical.toString(), null));
    }

    @Override
    public CategoryType category() { return GENE_COPY_NUMBER; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();

        // load only the genes in the driver catalog if the driver category is being run

        final Set<String> driverGenes = mConfig.DriverGenes;

        if(driverGenes.isEmpty())
            return items;

        try
        {
            List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(
                    fileSources.Purple, sampleId));

            copyNumbers.stream().filter(x -> driverGenes.contains(x.geneName())).forEach(x -> items.add(new GeneCopyNumberData(x, mFields)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read gene copy number data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }

    @Override
    public List<ComparableItem> loadFromExpected(final Map<String, List<ExpectedValue>> valuesByKey)
    {
        List<ComparableItem> items = Lists.newArrayList();

        for(List<ExpectedValue> expectedValues : valuesByKey.values())
        {
            GeneCopyNumberData geneCopyNumber = GeneCopyNumberData.fromExpected(expectedValues, mFields);
            geneCopyNumber.setRemark(ExpectedValue.remark(expectedValues));
            items.add(geneCopyNumber);
        }

        return items;
    }

    @Override
    public void compareItems(
            final List<Mismatch> mismatches, final MatchLevel matchLevel, final boolean includeMatches,
            final boolean includesExpected, final List<ComparableItem> oldItems, final List<ComparableItem> newItems)
    {
        if(mConfig.isExpectedSourced(SourceType.OLD))
            newItems.removeIf(x -> oldItems.stream().noneMatch(y -> y.matches(x)));

        CommonUtils.compareItems(this, mismatches, matchLevel, includeMatches, includesExpected, oldItems, newItems);
    }
}
