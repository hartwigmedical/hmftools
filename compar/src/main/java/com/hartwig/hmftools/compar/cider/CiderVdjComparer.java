package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.CDR3_SEQUENCE;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3Sequence;
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
        return List.of(Cdr3SequenceFile.Column.cdr3Seq.name(), Cdr3SequenceFile.Column.filter.name(), Cdr3SequenceFile.Column.locus.name());
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return new ArrayList<>();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = new ArrayList<>();

        try
        {
            List<Cdr3Sequence> cdr3Sequences = Cdr3SequenceFile.read(Cdr3SequenceFile.generateFilename(fileSources.Cider, sampleId));

            for(Cdr3Sequence cdr3Sequence : cdr3Sequences)
            {
                comparableItems.add(new CiderVdjData(cdr3Sequence));
            }
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load cider VDJ sequence data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
