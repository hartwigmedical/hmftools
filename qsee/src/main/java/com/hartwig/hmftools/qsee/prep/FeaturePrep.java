package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.qsee.cohort.FeatureMatrix;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class FeaturePrep
{
    private final CommonPrepConfig mConfig;

    public FeaturePrep(final CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SampleFeatures prepSample(SampleType sampleType, String sampleId)
    {
        QC_LOGGER.info("Extracting sample data - sampleType({}) sample({})", sampleType, sampleId);

        List<Feature> features = new ArrayList<>();

        List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mConfig).createCategoryPreps();
        for(CategoryPrep categoryPrep : categoryPreps)
        {
            QC_LOGGER.debug("sampleType({}) sample({}) - extracting category({})", sampleType, sampleId, categoryPrep.name());

            CategoryPrepTask task = new CategoryPrepTask(categoryPrep, sampleId, sampleType, mConfig.AllowMissingInput);

            task.run();

            List<Feature> categoryFeatures = task.getOutput();
            features.addAll(categoryFeatures);

        }

        return new SampleFeatures(sampleId, sampleType, features);
    }

    public FeatureMatrix prepCohortCategory(CategoryPrep categoryPrep, SampleType sampleType)
    {
        List<String> sampleIds = mConfig.getSampleIds(sampleType);

        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), sampleIds);

        List<Runnable> sampleCategoryTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < sampleIds.size(); ++sampleIndex)
        {
            CategoryPrepTask task = new CategoryPrepTask(
                    categoryPrep,
                    sampleIds.get(sampleIndex), sampleIndex, sampleIds.size(), sampleType,
                    sampleFeatureMatrix, mConfig.AllowMissingInput
            );

            sampleCategoryTasks.add(task);
        }

        TaskExecutor.executeRunnables(sampleCategoryTasks, mConfig.Threads);
        sampleCategoryTasks.clear();

        sampleFeatureMatrix.sortFeatureKeys();

        return sampleFeatureMatrix;
    }

    public List<SampleFeatures> prepCohort(SampleType sampleType)
    {
        List<SampleFeatures> cohortFeatures = new ArrayList<>();

        List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mConfig).createCategoryPreps();
        for(CategoryPrep categoryPrep : categoryPreps)
        {
            String logPrefix = CategoryPrepTask.logPrefix(sampleType, categoryPrep);

            QC_LOGGER.info("{} Extracting cohort data", logPrefix);

            FeatureMatrix sampleFeatureMatrix = prepCohortCategory(categoryPrep, sampleType);

            List<SampleFeatures> categorySampleFeatures = sampleFeatureMatrix.getRowIds().stream()
                    .map(sampleId -> getSampleFeatures(sampleFeatureMatrix, sampleId, sampleType))
                    .toList();

            cohortFeatures.addAll(categorySampleFeatures);
        }

        return cohortFeatures;
    }

    private static SampleFeatures getSampleFeatures(FeatureMatrix sampleFeatureMatrix, String sampleId, SampleType sampleType)
    {
        List<FeatureKey> featureKeys = sampleFeatureMatrix.getFeatureKeys();
        double[] featureValues = sampleFeatureMatrix.getRowValues(sampleId);

        List<Feature> features = IntStream.range(0, featureKeys.size())
                .mapToObj(i -> new Feature(featureKeys.get(i), featureValues[i]))
                .toList();

        return new SampleFeatures(sampleId, sampleType, features);
    }
}
