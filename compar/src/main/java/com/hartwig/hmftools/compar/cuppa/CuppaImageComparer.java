package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA_IMAGE;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CuppaImageComparer implements ItemComparer
{
    public static final String FLD_VIS_IMAGE = "cuppa_vis_image";

    private final ComparConfig mConfig;

    public CuppaImageComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return CUPPA_IMAGE; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_VIS_IMAGE, Double.NaN, 0);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = new ArrayList<>();

        String plotPath = CuppaPredictions.generateVisPlotFilename(fileSources.Cuppa, sampleId);
        comparableItems.add(new CuppaImageData(FLD_VIS_IMAGE, plotPath));

        return comparableItems;
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_VIS_IMAGE);
    }
}
