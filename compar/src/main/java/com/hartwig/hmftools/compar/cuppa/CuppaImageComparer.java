package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA_IMAGE;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ImageComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;

public class CuppaImageComparer extends ImageComparer
{
    public static final String IMAGE_NAME = "cuppa_vis_image";

    public CuppaImageComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config, 0.0, null, fieldCheckMap);
    }

    @Override
    public CategoryType category() { return CUPPA_IMAGE; }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = new ArrayList<>();
        String plotPath = CuppaPredictions.generateVisPlotFilename(fileSources.Cuppa, sampleId);
        CuppaImageData imageData = new CuppaImageData(IMAGE_NAME, plotPath);

        if(imageData.Image == null)
        {
            return null;
        }

        comparableItems.add(imageData);
        return comparableItems;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList();
    }
}
