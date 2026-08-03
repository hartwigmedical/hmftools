package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_LOCUS_SUMMARY;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummaryFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.IntField;

public class Cdr3LocusSummaryComparer extends ItemComparer
{
    protected static final String FLD_PASS_SEQUENCES = capitalize(Cdr3LocusSummaryFile.Column.passSequences.name());

    protected enum Fields
    {
        PassSequences;
    }

    public Cdr3LocusSummaryComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.PassSequences.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.PassSequences.toString(), null, 0.05),
                null));
    }

    @Override
    public CategoryType category()
    {
        return CDR3_LOCUS_SUMMARY;
    }

    /*
    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new IntField(FLD_PASS_SEQUENCES, i -> ((Cdr3LocusSummaryData) i).Cdr3LocusSummary.passSequences(),
                        true, null, 0.05)
        );
    }
    */

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(Fields.PassSequences.toString());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        try
        {
            return Cdr3LocusSummaryFile.read(Cdr3LocusSummaryFile.generateFilename(fileSources.Cider, sampleId))
                    .stream()
                    .map(x -> new Cdr3LocusSummaryData(x, mFields))
                    .collect(Collectors.toList());
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load cider locus summary data: {}", sampleId, e.toString());
            return null;
        }
    }
}
