package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.driver.DriverType.AMP;
import static com.hartwig.hmftools.common.driver.DriverType.DEL;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_PURPLE_SOMATIC_COPY_NUMBER;
import static com.hartwig.hmftools.common.driver.DriverType.HET_DEL;
import static com.hartwig.hmftools.common.driver.DriverType.PARTIAL_AMP;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
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
import com.hartwig.hmftools.compar.driver.DriverComparer;
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
}
