package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.GERMLINE_FLAGSTAT;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.determineFlagStatsFilePath;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SampleFileSources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineFlagstatComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineFlagstatComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return GERMLINE_FLAGSTAT;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_MAPPED_PROPORTION, MAPPED_PROPORTION_ABS_THRESHOLD, MAPPED_PROPORTION_PCT_THRESHOLD);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_MAPPED_PROPORTION);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final SampleFileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            BamFlagStats flagstat = BamFlagStats.read(determineFlagStatsFilePath(germlineSampleId, fileSources.germlineFlagstat()));
            comparableItems.add(new GermlineFlagstatData(flagstat));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load germline flagstat data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }
}
