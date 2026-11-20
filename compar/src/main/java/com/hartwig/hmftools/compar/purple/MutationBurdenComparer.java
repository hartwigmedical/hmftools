package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.MUTATION_BURDEN;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class MutationBurdenComparer extends ItemComparer
{
    protected enum Fields
    {
        TmbPerMb,
        MsIndelsPerMb,
        Tml,
        MsStatus,
        TmbStatus,
        TmlStatus;
    }

    public MutationBurdenComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.TmbPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TmbPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MsIndelsPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MsIndelsPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Tml.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Tml.toString(), 1.0, 0.05),
                null));

        mFields.add(new FieldInfo(
                Fields.MsStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.MsStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.TmbStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TmbStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.TmlStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TmlStatus.toString()),
                null));
    }

    @Override
    public CategoryType category() { return MUTATION_BURDEN; }

    @Override
    public boolean isMultiItem() { return false; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(
            final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            PurityContext purityContext = PurityContextFile.read(fileSources.Purple, sampleId);
            comparableItems.add(new MutationBurdenData(purityContext, mFields));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Purple mutation burden data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    @Override
    public List<ComparableItem> loadFromExpected(final Map<String,List<ExpectedValue>> valuesByKey)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        if(valuesByKey.size() == 1)
        {
            List<ExpectedValue> expectedValues = valuesByKey.get(MUTATION_BURDEN.toString());
            MutationBurdenData mutationBurden = new MutationBurdenData(expectedValues, mFields);
            mutationBurden.setRemark(ExpectedValue.remark(expectedValues));
            comparableItems.add(mutationBurden);
        }

        return comparableItems;
    }
}
