package com.hartwig.hmftools.compar.cobalt;

import static com.hartwig.hmftools.compar.common.Category.COBALT_RATIO;

import java.util.LinkedList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CobaltRatioComparer implements ItemComparer
{
    static final String REFERENCE_READ_COUNT = "ReferenceReadCount";
    static final String TUMOR_READ_COUNT = "TumorReadCount";
    static final String REFERENCE_GC_RATIO = "ReferenceGcRatio";
    static final String TUMOR_GC_RATIO = "TumorGcRatio";
    static final String REFERENCE_GC_DIPLOID_RATIO = "ReferenceGcDiploidRatio";
    private final ComparConfig mConfig;

    public CobaltRatioComparer(final ComparConfig mConfig)
    {
        this.mConfig = mConfig;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(REFERENCE_READ_COUNT, 0.0001, 0.0001);
        thresholds.addFieldThreshold(TUMOR_READ_COUNT, 0.0001, 0.0001);
        thresholds.addFieldThreshold(REFERENCE_GC_RATIO, 0.0001, 0.0001);
        thresholds.addFieldThreshold(TUMOR_GC_RATIO, 0.0001, 0.0001);
        thresholds.addFieldThreshold(REFERENCE_GC_DIPLOID_RATIO, 0.00001, 0.0001);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        String fileName = CobaltRatioFile.generateFilenameForReading(fileSources.Cobalt, sampleId);
        return new LinkedList<>(RawCobaltRatioFile.read(fileName));
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(REFERENCE_READ_COUNT, TUMOR_READ_COUNT, REFERENCE_GC_RATIO, TUMOR_GC_RATIO, REFERENCE_GC_DIPLOID_RATIO);
    }

    @Override
    public Category category()
    {
        return COBALT_RATIO;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        return List.of();
    }

    @Override
    public boolean hasReportable()
    {
        return false;
    }
}
