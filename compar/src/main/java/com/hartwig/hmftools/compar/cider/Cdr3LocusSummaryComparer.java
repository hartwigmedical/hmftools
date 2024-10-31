package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.CDR3_LOCUS_SUMMARY;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummaryFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class Cdr3LocusSummaryComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public Cdr3LocusSummaryComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return CDR3_LOCUS_SUMMARY; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(Cdr3LocusSummaryFile.Column.passSequences.name(), Double.NaN, 0.05);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Cdr3LocusSummaryData.comparedFieldNames();
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        return dbAccess.readCdr3LocusSummaries(sampleId).stream()
                .map(Cdr3LocusSummaryData::new)
                .collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        try
        {
            return Cdr3LocusSummaryFile.read(Cdr3LocusSummaryFile.generateFilename(fileSources.Cider, sampleId))
                    .stream()
                    .map(Cdr3LocusSummaryData::new)
                    .collect(Collectors.toList());
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load cider locus summary data: {}", sampleId, e.toString());
            return null;
        }
    }
}
