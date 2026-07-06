package com.hartwig.hmftools.compar;

import java.util.List;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public interface ItemComparer
{
    CategoryType category();

    boolean processSample(final String sampleId, final List<Mismatch> mismatches);

    List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType);

    List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources);

    List<Field> fields();

    List<String> displayFieldNames();

    default boolean hasReportable() { return true; }

}
