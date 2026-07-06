package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.cider.CiderVdjData.FILTER_FIELD;
import static com.hartwig.hmftools.compar.cider.CiderVdjData.LOCUS_FIELD;
import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CiderVdjComparer implements ItemComparer
{
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
    public List<Field> fields()
    {
        return List.of(
                new StringField(FILTER_FIELD, i -> ((CiderVdjData) i).mCdr3Sequence.filter(), true),
                new StringField(LOCUS_FIELD, i -> ((CiderVdjData) i).mCdr3Sequence.locus(), true)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(FILTER_FIELD, CiderVdjData.LOCUS_FIELD);
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
