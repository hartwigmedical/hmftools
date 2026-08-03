package com.hartwig.hmftools.compar.sigs;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class SigsComparer extends ItemComparer
{
    protected enum Fields
    {
        Percent;
    }

    public SigsComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Percent.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Percent.toString(), 0.05, null),
                "%.4f"));
    }

    @Override
    public CategoryType category()
    {
        return CategoryType.SIGS;
    }

    @Override
    public boolean hasReportable()
    {
        return false;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(Fields.Percent.toString());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            String filename = SignatureAllocationFile.generateFilename(fileSources.Sigs, sampleId);
            SignatureAllocationFile.read(filename).stream().map(x -> new SigsData(x, mFields)).forEach(comparableItems::add);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Sigs data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
