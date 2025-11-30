package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.V_CHORD;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_BREAST;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OTHER;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OVARIAN;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PANCREATIC;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PROSTATE;

import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.common.vchord.VChordPredictionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class VChordComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public VChordComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return V_CHORD; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_BREAST, 0.1, 0);
        thresholds.addFieldThreshold(FLD_OVARIAN, 0.1, 0);
        thresholds.addFieldThreshold(FLD_PANCREATIC, 0.1, 0);
        thresholds.addFieldThreshold(FLD_PROSTATE, 0.1, 0);
        thresholds.addFieldThreshold(FLD_OTHER, 0.1, 0);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_BREAST, FLD_OVARIAN, FLD_PANCREATIC, FLD_PROSTATE, FLD_OTHER);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Collections.emptyList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            VChordPrediction vChordData = VChordPredictionFile.read(VChordPredictionFile.generateFilename(fileSources.VChord, sampleId));
            comparableItems.add(new VChordData(vChordData));
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load vChord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
