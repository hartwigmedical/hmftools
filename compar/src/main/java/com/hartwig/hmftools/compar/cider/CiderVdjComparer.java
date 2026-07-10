package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
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
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CiderVdjComparer implements ItemComparer
{
    protected static final String FLD_FILTER = capitalize(Cdr3SequenceFile.Column.filter.name());
    protected static final String FLD_LOCUS = capitalize(Cdr3SequenceFile.Column.locus.name());

    private final ComparConfig mConfig;

    public CiderVdjComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() {
        return CDR3_SEQUENCE;
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_FILTER, i -> ((CiderVdjData) i).mCdr3Sequence.filter(), true),
                new StringField(FLD_LOCUS, i -> ((CiderVdjData) i).mCdr3Sequence.locus(), true)
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
        return List.of(FLD_FILTER, FLD_LOCUS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        return dbAccess.readCdr3Sequences(sampleId).stream()
                .map(CiderVdjData::new)
                .collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        try
        {
            return Cdr3SequenceFile.read(Cdr3SequenceFile.generateFilename(fileSources.Cider, sampleId)).stream()
                    .filter(seq -> seq.filter().equals("PASS") || seq.filter().equals("PARTIAL"))
                    .map(CiderVdjData::new)
                    .collect(Collectors.toList());
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load cider VDJ sequence data: {}", sampleId, e.toString());
            return null;
        }
    }
}
