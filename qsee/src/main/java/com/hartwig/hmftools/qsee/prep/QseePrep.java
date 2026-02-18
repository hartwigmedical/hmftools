package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.cohort.CohortPercentiles;
import com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile;
import com.hartwig.hmftools.qsee.cohort.FeaturePercentiles;
import com.hartwig.hmftools.qsee.cohort.PercentileTransformer;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class QseePrep
{
    private final QseePrepConfig mConfig;

    private static final String COL_FEATURE_VALUE = "FeatureValue";
    private static final String COL_PERCENTILE_IN_COHORT = "PctInCohort";
    private static final String COL_QC_STATUS = "QcStatus";
    private static final String COL_PLOT_LABEL = "PlotLabel";

    private static final String SAMPLE_ID_MULTI = "MULTI_SAMPLE";

    public QseePrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    private List<SampleFeatures> runFeaturePrepFor(SampleType sampleType)
    {
        List<String> sampleIds = mConfig.CommonPrep.getSampleIds(sampleType);

        boolean hasSampleType = !sampleIds.isEmpty();
        if(!hasSampleType)
        {
            return new ArrayList<>();
        }

        if(sampleIds.size() == 1)
        {
            SampleFeatures sampleFeatures = new FeaturePrep(mConfig.CommonPrep).prepSample(sampleType, sampleIds.get(0));
            return List.of(sampleFeatures);
        }
        else
        {
            return new FeaturePrep(mConfig.CommonPrep).prepMultiSample(sampleType);
        }
    }

    private List<VisSampleData> getVisSampleData(List<SampleFeatures> multiSampleFeatures, CohortPercentiles cohortPercentiles)
    {
        List<VisSampleData> visSampleData = new ArrayList<>();

        QC_LOGGER.info("Creating vis data entries");

        for(SampleFeatures sampleFeatures : multiSampleFeatures)
        {
            QC_LOGGER.debug("Creating vis data entries - sampleType({}) sample({})", sampleFeatures.sampleType(), sampleFeatures.sampleId());

            for(Feature feature : sampleFeatures.features())
            {
                 QC_LOGGER.trace("sampleType({}) sample({}) - transforming featureKey({}) featureValue({}) to percentile",
                         sampleFeatures.sampleId(), sampleFeatures.sampleType(), feature.key(), feature.value());

                FeaturePercentiles featurePercentiles = cohortPercentiles.getFeaturePercentiles(sampleFeatures.sampleType(), feature.key());

                double percentileInCohort = Double.NaN;
                if(featurePercentiles != null)
                {
                    PercentileTransformer transformer = featurePercentiles.transformer();
                    double featureValue = feature.value();
                    percentileInCohort = transformer.featureValueToPercentile(featureValue);
                }

                VisSampleData visData = new VisSampleData(
                        sampleFeatures.sampleId(),
                        sampleFeatures.sampleType(),
                        feature, percentileInCohort
                );

                visSampleData.add(visData);
            }
        }

        return visSampleData;
    }

    private void writeToFile(String outputFile, List<VisSampleData> visSampleDataEntries)
    {
        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            QC_LOGGER.info("Writing sample vis data to: {}", outputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(COL_SAMPLE_ID);
            header.add(COL_SAMPLE_TYPE);
            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);
            header.add(COL_FEATURE_VALUE);
            header.add(COL_PERCENTILE_IN_COHORT);
            header.add(COL_QC_STATUS);
            header.add(COL_PLOT_LABEL);

            writer.write(header.toString());
            writer.newLine();

            for(VisSampleData entry : visSampleDataEntries)
            {
                Feature feature = entry.feature();
                FeatureKey featureKey = feature.key();

                StringJoiner line = new StringJoiner(TSV_DELIM);
                line.add(entry.sampleId());
                line.add(entry.sampleType().name());
                line.add(featureKey.sourceTool().name());
                line.add(featureKey.type().name());
                line.add(featureKey.name());

                String featureValue = QseeFileCommon.DECIMAL_FORMAT.format(feature.value());
                line.add(featureValue);

                String percentileInCohort = QseeFileCommon.DECIMAL_FORMAT.format(entry.percentileInCohort());
                line.add(percentileInCohort);

                String qcStatus = "";
                if(feature.qcStatus() != null)
                {
                    qcStatus = feature.qcStatus().toString();
                }

                line.add(qcStatus);
                line.add(feature.plotLabel());

                writer.write(line.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write to file: {}", outputFile, e);
            System.exit(1);
        }
    }

    public static String formOutputFilename(String basePath, String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + "." + QSEE_FILE_ID + ".vis.features.tsv.gz";
    }

    public static String formOutputFilename(CommonPrepConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                SAMPLE_ID_MULTI;

        return formOutputFilename(config.OutputDir, sampleId);
    }

    public void run()
    {
        QC_LOGGER.info("Running {}", this.getClass().getSimpleName());

        CohortPercentiles cohortPercentiles = CohortPercentilesFile.read(mConfig.CohortPercentilesFile);

        List<SampleFeatures> multiSampleFeatures = new ArrayList<>();
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.TUMOR));
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.NORMAL));

        if(!mConfig.CommonPrep.isSinglePatient())
            multiSampleFeatures.sort(Comparator.comparing(SampleFeatures::sampleId));

        List<VisSampleData> visDataEntries = getVisSampleData(multiSampleFeatures, cohortPercentiles);
        String outputFile = formOutputFilename(mConfig.CommonPrep);
        writeToFile(outputFile, visDataEntries);
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
