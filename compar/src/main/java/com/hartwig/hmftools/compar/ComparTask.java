package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.driver.DriverType.AMP;
import static com.hartwig.hmftools.common.driver.DriverType.DEL;
import static com.hartwig.hmftools.common.driver.DriverType.PARTIAL_AMP;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_ERROR;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.InvalidDataItem;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.purple.GeneCopyNumberComparer;

public class ComparTask implements Callable<Void>
{
    private final int mTaskId;
    private final ComparConfig mConfig;
    private final List<String> mSampleIds;
    private final List<ItemComparer> mComparers;

    private final MismatchWriter mWriter;

    public ComparTask(int taskId, final ComparConfig config, final MismatchWriter writer)
    {
        mTaskId = taskId;
        mConfig = config;
        mWriter = writer;

        mSampleIds = Lists.newArrayList();
        mComparers = buildComparers(config);
    }

    public List<String> getSampleIds() { return mSampleIds; }

    @Override
    public Void call()
    {
        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            processSample(sampleId);

            if(i > 0 && (i % 50) == 0)
            {
                CMP_LOGGER.info("{}: processed {} samples", mTaskId, i);
            }
        }

        if(mConfig.Threads > 1)
        {
            CMP_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
        }

        return null;
    }

    private void processSample(final String sampleId)
    {
        int totalMismatches = 0;
        int failedTypes = 0;
        for(ItemComparer comparer : mComparers)
        {
            List<Mismatch> mismatches = Lists.newArrayList();

            try
            {
                if(mConfig.runCopyNumberGeneComparer() && comparer.category() == GENE_COPY_NUMBER)
                {
                    ((GeneCopyNumberComparer) comparer).addDriverGenes(loadCombinedCopyNumberDriverGenes(sampleId));
                }

                boolean status = comparer.processSample(sampleId, mismatches);

                if(!status)
                    ++failedTypes;
            }
            catch(Exception e)
            {
                CMP_LOGGER.error("sample({}) failed processing: {}", sampleId, e.toString());
                e.printStackTrace();
                ++failedTypes;

                InvalidDataItem invalidDataItem = new InvalidDataItem(comparer.category());
                mismatches = List.of(new Mismatch(invalidDataItem, null, INVALID_ERROR, Collections.emptyList()));
            }

            mWriter.writeSampleMismatches(sampleId, comparer, mismatches);
            totalMismatches += mismatches.size();
        }

        if(failedTypes == 0)
        {
            CMP_LOGGER.debug("sample({}) mismatches({})", sampleId, totalMismatches);
        }
        else
        {
            CMP_LOGGER.warn("sample({}) mismatches({}) failed types({})", sampleId, totalMismatches, failedTypes);
        }
    }

    private Set<String> loadCombinedCopyNumberDriverGenes(final String sampleId)
    {
        Set<String> combinedGenes = Sets.newHashSet();

        for(String sourceName : mConfig.SourceNames)
        {
            String sourceSampleId = mConfig.sourceSampleId(sourceName, sampleId);
            String sourceGermlineSampleId = mConfig.sourceGermlineSampleId(sourceName, sampleId);

            FileSources fileSources = FileSources.sampleInstance(mConfig.FileSources.get(sourceName), sourceSampleId, sourceGermlineSampleId);

            try
            {
                String purpleDriverFile = DriverCatalogFile.generateSomaticFilename(fileSources.Purple, sourceSampleId);

                DriverCatalogFile.read(purpleDriverFile).stream()
                        .filter(x -> x.driver() == AMP || x.driver() == PARTIAL_AMP || x.driver() == DEL)
                        .forEach(x -> combinedGenes.add(x.gene()));

            }
            catch(IOException e)
            {
                CMP_LOGGER.warn("sample({}) failed to load driver data: {}", sampleId, e.toString());
            }
        }

        return combinedGenes;
    }
}
