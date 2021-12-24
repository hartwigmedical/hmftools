package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CopyNumberComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public CopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return COPY_NUMBER; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<PurpleCopyNumber> copyNumbers = dbAccess.readCopynumbers(sampleId);

        List<ComparableItem> items = Lists.newArrayList();

        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            items.add(new CopyNumberData(copyNumber));
        }

        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(
                    PurpleCopyNumberFile.generateFilenameForReading(fileSources.Purple, sampleId));

            copyNumbers.forEach(x -> comparableItems.add(new CopyNumberData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Purple copy number data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }
}
