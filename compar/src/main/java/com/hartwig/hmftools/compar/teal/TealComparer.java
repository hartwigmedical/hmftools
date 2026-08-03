package com.hartwig.hmftools.compar.teal;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.TELOMERE_LENGTH;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.common.teal.TelomereLengthFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class TealComparer extends ItemComparer
{
    protected enum Fields
    {
        TelomereLength;
    }

    public TealComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.TelomereLength.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TelomereLength.toString(), null, 0.05),
                "%.2f"));
    }

    @Override
    public CategoryType category() { return TELOMERE_LENGTH; }

    // TODO: remove
//    @Override
//    public List<Field> fields(final MatchLevel matchLevel)
//    {
//        return List.of(
//                new DoubleField(FLD_TELOMERE_LENGTH, i -> ((TealData) i).TelomereLength.finalTelomereLength(),
//                        true, null, 0.05, "%.2f")
//        );
//    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(Fields.TelomereLength.toString());
    }


    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        try
        {
            TelomereLength telomereLength = TelomereLengthFile.read(TelomereLengthFile.generateFilename(fileSources.Teal, sampleId));
            return Lists.newArrayList(new TealData(telomereLength, mFields));
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Teal data: {}", sampleId, e.toString());
            return null;
        }
    }
}
