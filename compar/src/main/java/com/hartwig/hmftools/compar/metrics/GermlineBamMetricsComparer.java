package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.GERMLINE_BAM_METRICS;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.metrics.GermlineBamMetricsData.FLD_PERCENTAGE_10X;
import static com.hartwig.hmftools.compar.metrics.GermlineBamMetricsData.FLD_PERCENTAGE_20X;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.DUPLICATE_PERCENTAGE_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_DUPLICATE_PERCENTAGE;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.OLD_BAM_METRICS_FILE_EXTENSION;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineBamMetricsComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineBamMetricsComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return GERMLINE_BAM_METRICS;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_DUPLICATE_PERCENTAGE, DUPLICATE_PERCENTAGE_ABS_THRESHOLD, DUPLICATE_PERCENTAGE_PCT_THRESHOLD);
        thresholds.addFieldThreshold(FLD_PERCENTAGE_10X, 0.03, 0);
        thresholds.addFieldThreshold(FLD_PERCENTAGE_20X, 0.03, 0);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_DUPLICATE_PERCENTAGE, FLD_PERCENTAGE_10X, FLD_PERCENTAGE_20X);
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
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            WGSMetrics metrics = WGSMetricsFile.read(determineFilePath(germlineSampleId, fileSources));
            comparableItems.add(new GermlineBamMetricsData(metrics));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load germline BAM metrics data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

    private static String determineFilePath(final String germlineSampleId, final FileSources fileSources)
    {
        String currentFileName = WGSMetricsFile.generateFilename(fileSources.GermlineBamMetrics, germlineSampleId);
        String oldFileName = checkAddDirSeparator(fileSources.GermlineBamMetrics) + germlineSampleId + OLD_BAM_METRICS_FILE_EXTENSION;
        if(!fileExists(currentFileName) && fileExists(oldFileName))
        {
            return oldFileName;
        }
        else
        {
            return currentFileName;
        }
    }
}
