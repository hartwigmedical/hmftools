package com.hartwig.hmftools.compar.cuppa;

import static com.hartwig.hmftools.common.cuppa.DataType.PROB;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.CUPPA;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_PROBABILITY;
import static com.hartwig.hmftools.compar.cuppa.CuppaData.FLD_TOP_CANCER_TYPE;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.cuppa.ClassifierGroup;
import com.hartwig.hmftools.common.cuppa.ClassifierName;
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

import org.apache.logging.log4j.util.Strings;

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
    public List<ComparableItem> loadFromOrangeJson(final String sampleId, final JsonObject json)
    {
        JsonObject bestPredictionJson = json.getAsJsonObject("cuppa").getAsJsonObject("bestPrediction");
        String cancerType = bestPredictionJson.get("cancerType").getAsString();

        final List<ComparableItem> comparableItems = new ArrayList<>();
        comparableItems.add(createCuppaData(sampleId, ClassifierGroup.DNA, ClassifierName.GEN_POS, cancerType,
                bestPredictionJson.get("genomicPositionClassifier").getAsDouble(), 3));
        comparableItems.add(createCuppaData(sampleId, ClassifierGroup.DNA, ClassifierName.SNV96, cancerType,
                bestPredictionJson.get("snvPairwiseClassifier").getAsDouble(), 4));
        comparableItems.add(createCuppaData(sampleId, ClassifierGroup.DNA, ClassifierName.EVENT, cancerType,
                bestPredictionJson.get("featureClassifier").getAsDouble(), 5));
        if(bestPredictionJson.get("expressionPairwiseClassifier").isJsonNull())
        {
            comparableItems.add(createCuppaData(sampleId, ClassifierGroup.RNA, ClassifierName.DNA_COMBINED,  cancerType,
                    bestPredictionJson.get("likelihood").getAsDouble(), 1));
        }
        else
        {
            comparableItems.add(createCuppaData(sampleId, ClassifierGroup.RNA, ClassifierName.COMBINED, cancerType,
                    bestPredictionJson.get("likelihood").getAsDouble(), 1));
            comparableItems.add(createCuppaData(sampleId, ClassifierGroup.RNA, ClassifierName.GENE_EXP, cancerType,
                    bestPredictionJson.get("expressionPairwiseClassifier").getAsDouble(), 6));
            comparableItems.add(createCuppaData(sampleId, ClassifierGroup.RNA, ClassifierName.ALT_SJ, cancerType,
                    bestPredictionJson.get("altSjCohortClassifier").getAsDouble(), 7));
        }
        return comparableItems;
    }

    private static CuppaData createCuppaData(final String sampleId, final ClassifierGroup classifierGroup,
            final ClassifierName classifierName, final String cancerType, final double value, final int rankGroup)
    {
        final CuppaPredictionEntry predictionEntry = new CuppaPredictionEntry(sampleId, PROB, classifierGroup, classifierName,
                Strings.EMPTY, Double.NaN, cancerType, value, 1, rankGroup);
        return new CuppaData(predictionEntry);
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
