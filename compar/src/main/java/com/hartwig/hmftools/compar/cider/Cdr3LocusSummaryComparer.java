package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_LOCUS_SUMMARY;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummaryFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class Cdr3LocusSummaryComparer implements ItemComparer
{
    protected static final String FLD_PASS_SEQUENCES = capitalize(Cdr3LocusSummaryFile.Column.passSequences.name());

    private final ComparConfig mConfig;

    public Cdr3LocusSummaryComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category()
    {
        return CDR3_LOCUS_SUMMARY;
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new IntField(FLD_PASS_SEQUENCES, i -> ((Cdr3LocusSummaryData) i).Cdr3LocusSummary.passSequences(),
                        true, null, 0.05)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(FLD_PASS_SEQUENCES);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
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
