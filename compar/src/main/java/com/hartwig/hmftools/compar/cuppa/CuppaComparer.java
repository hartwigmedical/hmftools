package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.common.cuppa.DataType.PROB;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class CuppaComparer extends ItemComparer
{
    protected enum Fields
    {
        TopCancerType,
        Probability;
    }

    public CuppaComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.TopCancerType.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TopCancerType.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Probability.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Probability.toString(), 0.1, null),
                "%.3f"));
    }

    @Override
    public CategoryType category() { return CUPPA; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = new ArrayList<>();

        try
        {
            String visDataPath = CuppaPredictions.generateVisDataTsvFilename(fileSources.Cuppa, sampleId);

            CuppaPredictions topProbabilities = CuppaPredictions
                    .fromTsv(visDataPath)
                    .subsetByDataType(PROB)
                    .getTopPredictions(1);

            for(CuppaPredictionEntry predictionEntry : topProbabilities.PredictionEntries)
            {
                comparableItems.add(new CuppaData(predictionEntry, mFields));
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Cuppa data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
