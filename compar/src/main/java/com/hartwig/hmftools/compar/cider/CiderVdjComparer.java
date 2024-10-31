package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.CDR3_SEQUENCE;

import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CiderVdjComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public CiderVdjComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return CDR3_SEQUENCE; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return CiderVdjData.comparedFieldNames();
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        return dbAccess.readCdr3Sequences(sampleId).stream()
                .map(CiderVdjData::new)
                .collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
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
