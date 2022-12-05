package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.common.cuppa.ClassifierType.GENDER;
import static com.hartwig.hmftools.common.cuppa.CuppaDataFile.getRankedCancerTypes;
import static com.hartwig.hmftools.compar.Category.CUPPA;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_TOP_CANCER_TYPE;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.ClassifierType;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.ResultType;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.driver.DriverData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CuppaComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public CuppaComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return CUPPA; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_LIKELIHOOD, 0.1, 0);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_TOP_CANCER_TYPE, FLD_LIKELIHOOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<CuppaDataFile> cuppaDataList = CuppaDataFile.read(CuppaDataFile.generateFilename(fileSources.Cuppa, sampleId));

            for(CuppaDataFile cuppaData : cuppaDataList)
            {
                if(cuppaData.Result != ResultType.CLASSIFIER || cuppaData.DataType.equals(GENDER.toString()))
                    continue;

                if(cuppaData.CancerTypeValues.isEmpty())
                    continue;

                List<String> rankedCancerTypes = getRankedCancerTypes(cuppaData.CancerTypeValues);
                String topRefCancerType = rankedCancerTypes.get(0);
                double topCancerValue = cuppaData.CancerTypeValues.get(topRefCancerType);
                comparableItems.add(new CuppaData(new ClassifierData(cuppaData.DataType, topRefCancerType, topCancerValue)));
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Cuppa data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
