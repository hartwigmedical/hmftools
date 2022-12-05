package com.hartwig.hmftools.compar;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public interface ItemComparer
{
    Category category();

    boolean processSample(final String sampleId, final List<Mismatch> mismatches);

    List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess);

    List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources);

    List<String> comparedFieldNames();

    default void registerThresholds(final DiffThresholds thresholds) {}

    default boolean hasReportable() { return true; }

}
