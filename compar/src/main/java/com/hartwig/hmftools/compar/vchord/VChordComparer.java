package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.V_CHORD;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.common.vchord.VChordPredictionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class VChordComparer extends ItemComparer
{
    protected enum Fields
    {
        BreastCancerHrdScore,
        OvarianCancerHrdScore,
        OtherCancerScore;
    }

    public VChordComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.BreastCancerHrdScore.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.BreastCancerHrdScore.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.OvarianCancerHrdScore.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.OvarianCancerHrdScore.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.OtherCancerScore.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.OtherCancerScore.toString(), 0.1, null),
                "%.2f"));
    }

    @Override
    public CategoryType category() { return V_CHORD; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            VChordPrediction vChordData = VChordPredictionFile.read(VChordPredictionFile.generateFilename(fileSources.VChord, sampleId));
            comparableItems.add(new VChordData(vChordData, mFields));
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load vChord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
