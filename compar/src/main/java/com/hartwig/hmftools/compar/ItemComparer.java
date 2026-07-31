package com.hartwig.hmftools.compar;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.FieldCheck;

public interface ItemComparer
{
    CategoryType category();

    boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldCheckCache fieldConfig);

    default List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        return Collections.emptyList();
    }

    default List<ComparableItem> loadFromFile(
            final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources,
            final Map<String,FieldCheck> fieldCheckMap)
    {
        return loadFromFile(sampleId, germlineSampleId, fileSources);
    }

    default List<ComparableItem> loadFromTruthset(final List<TruthsetValue> truthsetValues, final Map<String,FieldCheck> fieldCheckMap)
    {
        return Collections.emptyList();
    }

    List<Field> fields(final MatchLevel matchLevel);

    List<String> displayFieldNames();

    default boolean hasReportable() { return true; }

}
