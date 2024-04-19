package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_MAJOR_ALLELE_CN;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_METHOD;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
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
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_COPY_NUMBER, 0.5, 0.15);
        thresholds.addFieldThreshold(FLD_MAJOR_ALLELE_CN, 0.5, 0.15);
    }

    @Override
    public boolean hasReportable() { return false; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_COPY_NUMBER, FLD_MAJOR_ALLELE_CN, FLD_METHOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<PurpleCopyNumber> copyNumbers = dbAccess.readCopynumbers(sampleId);
        List<ComparableItem> items = Lists.newArrayList();
        copyNumbers.forEach(x -> items.add(new CopyNumberData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(
                    fileSources.Purple, sampleId));

            copyNumbers.forEach(x -> comparableItems.add(new CopyNumberData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read copy number data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
