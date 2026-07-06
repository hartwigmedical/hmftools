package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA_IMAGE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ImageComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CuppaImageComparer extends ImageComparer
{
    public static final String FLD_VIS_IMAGE = "cuppa_vis_image";

    private final ComparConfig mConfig;

    public CuppaImageComparer(final ComparConfig config)
    {
        super(null, 0.);
        mConfig = config;
    }

    @Override
    public CategoryType category() { return CUPPA_IMAGE; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        String plotPath = CuppaPredictions.generateVisPlotFilename(fileSources.Cuppa, sampleId);
        CuppaImageData imageData = new CuppaImageData(FLD_VIS_IMAGE, plotPath);
        if(imageData.Image != null)
        {
            return List.of(imageData);
        }
        else
        {
            return null;
        }
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(FLD_VIS_IMAGE);
    }
}
