package com.hartwig.hmftools.compar;

import java.util.List;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.Field;

public interface ItemComparer
{
    CategoryType category();

    boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig);

    List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources);

    List<Field> fields(final MatchLevel matchLevel);

    List<String> displayFieldNames();

    default boolean hasReportable() { return true; }

}
