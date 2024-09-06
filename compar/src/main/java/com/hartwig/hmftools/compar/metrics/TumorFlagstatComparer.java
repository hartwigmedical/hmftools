package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.TUMOR_FLAGSTAT;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.FLD_MAPPED_PROPORTION;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.MAPPED_PROPORTION_PCT_THRESHOLD;
import static com.hartwig.hmftools.compar.metrics.MetricsCommon.OLD_FLAGSTAT_FILE_EXTENSION;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class TumorFlagstatComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public TumorFlagstatComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return TUMOR_FLAGSTAT;
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
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            Flagstat flagstat = FlagstatFile.read(determineFilePath(sampleId, fileSources));
            comparableItems.add(new TumorFlagstatData(flagstat));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load tumor flagstat data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

    private static String determineFilePath(final String sampleId, final FileSources fileSources)
    {
        String currentFileName = FlagstatFile.generateFilename(fileSources.TumorFlagstat, sampleId);
        String oldFileName = checkAddDirSeparator(fileSources.TumorFlagstat) + sampleId + OLD_FLAGSTAT_FILE_EXTENSION;
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
