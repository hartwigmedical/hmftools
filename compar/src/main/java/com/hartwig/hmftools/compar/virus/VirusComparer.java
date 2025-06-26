package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.VIRUS_INTERPRETER_DIR;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.VIRUS;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_DRIVER_LIKELIHOOD;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_INTEGRATIONS;
import static com.hartwig.hmftools.compar.virus.VirusData.FLD_MEAN_COVERAGE;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SampleFileSources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;

public class VirusComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public VirusComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return VIRUS;
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_MEAN_COVERAGE, 0, 0.15);
        thresholds.addFieldThreshold(FLD_INTEGRATIONS, 0, 0.20);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_REPORTED, FLD_INTEGRATIONS, FLD_MEAN_COVERAGE, FLD_DRIVER_LIKELIHOOD);
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
            AnnotatedVirusFile.read(determineFileName(sampleId, fileSources))
                    .forEach(v -> comparableItems.add(new VirusData(v)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Virus interpreter data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

    @NotNull
    private static String determineFileName(final String sampleId, final SampleFileSources fileSources)
    {
        // dirty hack to get old virus directory working automatically most of the time
        final String currentFileName = AnnotatedVirusFile.generateFileName(fileSources.virus(), sampleId);
        final String oldFileName =
                AnnotatedVirusFile.generateFileName(fileSources.virus().replaceAll(VIRUS_INTERPRETER_DIR, "virus_interpreter"), sampleId);

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
