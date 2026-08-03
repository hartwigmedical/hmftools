package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class CiderVdjComparer extends ItemComparer
{
    protected enum Fields
    {
        Filter,
        Locus;
    }

    public CiderVdjComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Filter.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Filter.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Locus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Locus.toString()), null));
    }

    @Override
    public CategoryType category() {
        return CDR3_SEQUENCE;
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.name()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        try
        {
            return Cdr3SequenceFile.read(Cdr3SequenceFile.generateFilename(fileSources.Cider, sampleId)).stream()
                    .filter(seq -> seq.filter().equals("PASS") || seq.filter().equals("PARTIAL"))
                    .map(seq -> new CiderVdjData(seq, mFields))
                    .collect(Collectors.toList());
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load cider VDJ sequence data: {}", sampleId, e.toString());
            return null;
        }
    }
}
