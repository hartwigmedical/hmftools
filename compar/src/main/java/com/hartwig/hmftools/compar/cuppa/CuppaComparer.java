package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.common.cuppa.DataType.PROB;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.CUPPA;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_CLASSIFIER_NAME;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_PROBABILITY;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_TOP_CANCER_TYPE;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
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
        thresholds.addFieldThreshold(FLD_PROBABILITY, 0.1, 0);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_TOP_CANCER_TYPE, FLD_PROBABILITY);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = new ArrayList<>();

        try
        {
            String visDataPath = CuppaPredictions.generateVisDataTsvFilename(fileSources.Cuppa, sampleId);

            CuppaPredictions topProbabilities = CuppaPredictions
                    .fromTsv(visDataPath)
                    .subsetByDataType(PROB)
                    .getTopPredictions(1);

            for(CuppaPredictionEntry predictionEntry : topProbabilities.PredictionEntries)
            {
                comparableItems.add(new CuppaData(predictionEntry));
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
