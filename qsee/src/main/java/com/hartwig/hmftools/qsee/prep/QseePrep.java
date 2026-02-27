package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.cohort.CohortPercentiles;
import com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;

import org.jetbrains.annotations.Nullable;

public class QseePrep
{
    private final QseePrepConfig mConfig;

    public QseePrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    private List<SampleFeatures> runFeaturePrepFor(SampleType sampleType)
    {
        List<String> sampleIds = mConfig.getSampleIds(sampleType);

        boolean hasSampleType = !sampleIds.isEmpty();
        if(!hasSampleType)
        {
            return List.of();
        }

        if(sampleIds.size() == 1)
        {
            SampleFeatures sampleFeatures = new FeaturePrep(mConfig).prepSample(sampleType, sampleIds.get(0));
            return List.of(sampleFeatures);
        }
        else
        {
            return new FeaturePrep(mConfig).prepMultiSample(sampleType);
        }
    }

    private List<VisSampleData> getVisSampleData(List<SampleFeatures> multiSampleFeatures, @Nullable CohortPercentiles cohortPercentiles)
    {
        List<VisSampleData> visSampleData = new ArrayList<>();

        for(SampleFeatures sampleFeatures : multiSampleFeatures)
        {
            QC_LOGGER.info("Creating vis data entries - sampleType({}) sample({})",
                    sampleFeatures.sampleType(), sampleFeatures.sampleId());

            for(Feature feature : sampleFeatures.features())
            {
                if(cohortPercentiles != null)
                    cohortPercentiles.warnIfMissing(sampleFeatures.sampleType(), feature.key());

                VisSampleData visData = new VisSampleData(sampleFeatures.sampleId(), sampleFeatures.sampleType(), feature);
                visSampleData.add(visData);
            }
        }

        return visSampleData;
    }

    public void run()
    {
        QC_LOGGER.info("Running {}", this.getClass().getSimpleName());

        List<SampleFeatures> multiSampleFeatures = new ArrayList<>();
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.TUMOR));
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.NORMAL));

        if(!mConfig.isSinglePatient())
            multiSampleFeatures.sort(Comparator.comparing(SampleFeatures::sampleId));

        CohortPercentiles cohortPercentiles = (mConfig.CohortPercentilesFile != null)
                ? CohortPercentilesFile.read(mConfig.CohortPercentilesFile)
                : null;

        List<VisSampleData> visDataEntries = getVisSampleData(multiSampleFeatures, cohortPercentiles);

        String outputFile = VisDataFile.generateFilename(mConfig);
        VisDataFile.write(outputFile, visDataEntries);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseePrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        QseePrepConfig config = new QseePrepConfig(configBuilder);
        new QseePrep(config).run();
    }
}
