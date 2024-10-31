package com.hartwig.hmftools.compar.teal;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.TELOMERE_LENGTH;
import static com.hartwig.hmftools.compar.teal.TealData.FLD_TELOMERE_LENGTH;

import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.common.teal.TelomereLengthFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class TealComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public TealComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return TELOMERE_LENGTH; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_TELOMERE_LENGTH, Double.NaN, 0.05);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return TealData.comparedFieldNames();
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Collections.emptyList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            TelomereLength telomereLength = TelomereLengthFile.read(TelomereLengthFile.generateFilename(fileSources.Teal, sampleId));
            return Lists.newArrayList(new TealData(telomereLength));
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Teal data: {}", sampleId, e.toString());
            return null;
        }
    }
}
