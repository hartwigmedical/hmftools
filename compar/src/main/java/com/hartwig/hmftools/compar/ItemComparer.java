package com.hartwig.hmftools.compar;

import java.util.List;

import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public interface ItemComparer
{
    Category category();

    boolean processSample(final String sampleId, final List<Mismatch> mismatches);

    List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName);

    List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources);

    List<String> comparedFieldNames();

    default void registerThresholds(final DiffThresholds thresholds) {}

    default boolean hasReportable() { return true; }

}
