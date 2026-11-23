package com.hartwig.hmftools.qsee;

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
import com.hartwig.hmftools.qsee.feature.FeaturePrep;
import com.hartwig.hmftools.qsee.feature.SampleFeatures;
import com.hartwig.hmftools.qsee.lod.GeneLodCalculator;
import com.hartwig.hmftools.qsee.plot.QseePlot;
import com.hartwig.hmftools.qsee.plot.QseePlotConfig;
import com.hartwig.hmftools.qsee.plot.VisDataFile;
import com.hartwig.hmftools.qsee.status.QcStatusFile;

public class QseeApplication
{
    private final ConfigBuilder mConfigBuilder;
    private final QseeConfig mConfig;

    public QseeApplication(ConfigBuilder configBuilder)
    {
        mConfigBuilder = configBuilder;
        mConfig = new QseeConfig(configBuilder);
    }

    private List<SampleFeatures> prepSampleFeatures(SampleType sampleType)
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

    private static void warnIfMissingCohortPercentiles(List<SampleFeatures> multiSampleFeatures, CohortPercentiles cohortPercentiles)
    {
        for(SampleFeatures sampleFeatures : multiSampleFeatures)
        {
            for(Feature feature : sampleFeatures.features())
            {
                cohortPercentiles.warnIfMissing(sampleFeatures.sampleType(), feature.key());
            }
        }
    }

    public void run()
    {
        QC_LOGGER.info("Running {}", APP_NAME);

        List<SampleFeatures> multiSampleFeatures = new ArrayList<>();
        multiSampleFeatures.addAll(prepSampleFeatures(SampleType.TUMOR));
        multiSampleFeatures.addAll(prepSampleFeatures(SampleType.NORMAL));

        if(!mConfig.isSinglePatient())
            multiSampleFeatures.sort(Comparator.comparing(SampleFeatures::sampleId));

        if(mConfig.CohortPercentilesFile != null)
        {
            CohortPercentiles cohortPercentiles = CohortPercentilesFile.read(mConfig.CohortPercentilesFile);
            warnIfMissingCohortPercentiles(multiSampleFeatures, cohortPercentiles);
        }

        String visDataFile = VisDataFile.generateFilename(mConfig);
        VisDataFile.write(visDataFile, multiSampleFeatures);

        String qcStatusFile = QcStatusFile.generateFilename(mConfig);
        QcStatusFile.write(qcStatusFile, multiSampleFeatures);

        GeneLodCalculator geneLodCalculator = new GeneLodCalculator(mConfig);
        geneLodCalculator.run();

        if(mConfig.SkipPlot)
        {
            QC_LOGGER.info("Skipping plotting");
            return;
        }

        QseePlotConfig plotConfig = new QseePlotConfig(mConfigBuilder);
        QseePlot qseePlot = new QseePlot(plotConfig);
        qseePlot.run();
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseeConfig.registerConfig(configBuilder);
        QseePlotConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new QseeApplication(configBuilder).run();
    }
}
