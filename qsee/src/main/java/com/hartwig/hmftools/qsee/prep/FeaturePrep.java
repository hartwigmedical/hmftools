package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.qsee.cohort.FeatureMatrix;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.prep.category.BaseQualRecalibrationPrep;
import com.hartwig.hmftools.qsee.prep.category.CoverageDistributionPrep;
import com.hartwig.hmftools.qsee.prep.category.DiscordantFragFreqPrep;
import com.hartwig.hmftools.qsee.prep.category.DuplicateFreqPrep;
import com.hartwig.hmftools.qsee.prep.category.FragLengthDistributionPrep;
import com.hartwig.hmftools.qsee.prep.category.GcBiasPrep;
import com.hartwig.hmftools.qsee.prep.category.MissedGeneVariantPrep;
import com.hartwig.hmftools.qsee.prep.category.MsIndelErrorPrep;
import com.hartwig.hmftools.qsee.prep.category.SummaryTablePrep;

public class FeaturePrep
{
    private final CommonPrepConfig mConfig;

    public FeaturePrep(final CommonPrepConfig config)
    {
        mConfig = config;
    }

    public static List<CategoryPrep> createCategoryPreps(CommonPrepConfig config)
    {
        return List.of(
                new SummaryTablePrep(config),
                new CoverageDistributionPrep(config),
                new FragLengthDistributionPrep(config),
                new MissedGeneVariantPrep(config),
                new DuplicateFreqPrep(config),
                new GcBiasPrep(config),
                new DiscordantFragFreqPrep(config),
                new BaseQualRecalibrationPrep(config),
                new MsIndelErrorPrep(config)
        );
    }

    public SampleFeatures prepSample(SampleType sampleType, String sampleId)
    {
        QC_LOGGER.info("Extracting sample data - sampleType({}) sample({})", sampleType, sampleId);

        List<Feature> features = new ArrayList<>();

        List<CategoryPrep> categoryPreps = createCategoryPreps(mConfig);
        for(CategoryPrep categoryPrep : categoryPreps)
        {
            QC_LOGGER.debug("Extracting category({})", categoryPrep.name());

            CategoryPrepTask task = new CategoryPrepTask(categoryPrep, sampleId, sampleType, mConfig.AllowMissingInput);

            task.run();

            List<Feature> categoryFeatures = task.getOutput();
            features.addAll(categoryFeatures);

        }

        return new SampleFeatures(sampleId, sampleType, features);
    }

    public List<SampleFeatures> prepMultiSample(SampleType sampleType)
    {
        QC_LOGGER.info("Extracting multi-sample data - sampleType({})", sampleType.name());

        List<String> sampleIds = mConfig.getSampleIds(sampleType);

        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new HashMap<>(), mConfig.getSampleIds(sampleType));
        List<CategoryPrep> categoryPreps = createCategoryPreps(mConfig);
        for(CategoryPrep categoryPrep : categoryPreps)
        {
            QC_LOGGER.info("Extracting category({})", categoryPrep.name());
            prepCohortCategory(categoryPrep, sampleType, sampleFeatureMatrix, false);
        }

        sampleFeatureMatrix.sortFeatureKeys();

        List<SampleFeatures> multiSampleFeatures = new ArrayList<>();
        for(String sampleId : sampleIds)
        {
            List<FeatureKey> featureKeys = sampleFeatureMatrix.getFeatureKeys();
            double[] featureValues = sampleFeatureMatrix.getRowValues(sampleId);

            List<Feature> features = IntStream.range(0, featureKeys.size())
                    .mapToObj(featureIndex -> new Feature(featureKeys.get(featureIndex), featureValues[featureIndex]))
                    .toList();

            SampleFeatures sampleFeatures = new SampleFeatures(sampleId, sampleType, features);
            multiSampleFeatures.add(sampleFeatures);
        }

        return multiSampleFeatures;
    }

    public void prepCohortCategory(CategoryPrep categoryPrep, SampleType sampleType, FeatureMatrix outputMatrix, boolean isTraining)
    {
        List<String> sampleIds = mConfig.getSampleIds(sampleType);

        List<Runnable> sampleCategoryTasks = new ArrayList<>();
        AtomicInteger samplesMissingInputCount = new AtomicInteger(0);

        for(int sampleIndex = 0; sampleIndex < sampleIds.size(); ++sampleIndex)
        {
            CategoryPrepTask task = new CategoryPrepTask(
                    categoryPrep,
                    sampleIds.get(sampleIndex), sampleIndex, sampleIds.size(), sampleType, outputMatrix,
                    mConfig.AllowMissingInput, samplesMissingInputCount
            );

            sampleCategoryTasks.add(task);
        }

        TaskExecutor.executeRunnables(sampleCategoryTasks, mConfig.Threads);
        sampleCategoryTasks.clear();

        if(isTraining)
        {
            if(samplesMissingInputCount.get() == sampleIds.size())
            {
                QC_LOGGER.error("failed prep as no samples had data for sampleType({}) category({})",
                        sampleType, categoryPrep.name());

                System.exit(1);
            }
            else if(samplesMissingInputCount.get() > 0)
            {
                QC_LOGGER.warn("sampleType({}) category({}) - {}/{} samples had missing input files",
                        sampleType, categoryPrep.name(), samplesMissingInputCount.get(), sampleIds.size());
            }
        }
    }
}
