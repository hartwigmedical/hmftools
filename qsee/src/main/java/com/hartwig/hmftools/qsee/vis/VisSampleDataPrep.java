package com.hartwig.hmftools.qsee.vis;

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
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.prep.FeaturePrep;
import com.hartwig.hmftools.qsee.prep.SampleFeatures;

public class VisSampleDataPrep
{
    private final VisPrepConfig mVisPrepConfig;

    private static final String COL_PERCENTILE_IN_COHORT = "PercentileInCohort";
    private static final String SAMPLE_ID_MULTI = "MULTI_SAMPLE";

    public VisSampleDataPrep(VisPrepConfig config)
    {
        mVisPrepConfig = config;
    }

    private List<SampleFeatures> runFeaturePrepFor(SampleType sampleType)
    {
        List<String> sampleIds = mVisPrepConfig.CommonPrep.getSampleIds(sampleType);

        boolean hasSampleType = !sampleIds.isEmpty();
        if(!hasSampleType)
        {
            return new ArrayList<>();
        }

        if(sampleIds.size() == 1)
        {
            SampleFeatures sampleFeatures = new FeaturePrep(mVisPrepConfig.CommonPrep).prepSample(sampleType, sampleIds.get(0));
            return List.of(sampleFeatures);
        }
        else
        {
            return new FeaturePrep(mVisPrepConfig.CommonPrep).prepMultiSample(sampleType);
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
                PercentileTransformer transformer = featurePercentiles.transformer();
                double percentileInCohort = transformer.featureValueToPercentile(feature.value());

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
            header.add(COL_PERCENTILE_IN_COHORT);

            writer.write(header.toString());
            writer.newLine();

            for(VisSampleData entry : visSampleDataEntries)
            {
                FeatureKey featureKey = entry.feature().key();

                StringJoiner line = new StringJoiner(TSV_DELIM);
                line.add(entry.sampleId());
                line.add(entry.sampleType().name());
                line.add(featureKey.sourceTool().name());
                line.add(featureKey.type().name());
                line.add(featureKey.name());

                String percentileInCohort = CohortPercentilesFile.PERCENTILE_FORMAT.format(entry.percentileInCohort());
                line.add(percentileInCohort);

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

    private boolean isSinglePatient()
    {
        return mVisPrepConfig.CommonPrep.getSampleIds(SampleType.TUMOR).size() <= 1 &&
                mVisPrepConfig.CommonPrep.getSampleIds(SampleType.NORMAL).size() <= 1;
    }

    public void run()
    {
        QC_LOGGER.info("Running {}", this.getClass().getSimpleName());

        CohortPercentiles cohortPercentiles = CohortPercentilesFile.read(mVisPrepConfig.CohortPercentilesFile);

        List<SampleFeatures> multiSampleFeatures = new ArrayList<>();
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.TUMOR));
        multiSampleFeatures.addAll(runFeaturePrepFor(SampleType.NORMAL));

        if(!isSinglePatient())
            multiSampleFeatures.sort(Comparator.comparing(SampleFeatures::sampleId));

        List<VisSampleData> visDataEntries = getVisSampleData(multiSampleFeatures, cohortPercentiles);
        String sampleId = isSinglePatient() ? mVisPrepConfig.CommonPrep.getSampleIds(SampleType.TUMOR).get(0) : SAMPLE_ID_MULTI;
        String outputFile = checkAddDirSeparator(mVisPrepConfig.CommonPrep.OutputDir) + sampleId + "." + QSEE_FILE_ID + ".vis.features.tsv.gz";
        writeToFile(outputFile, visDataEntries);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        VisPrepConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        VisPrepConfig config = new VisPrepConfig(configBuilder);

        new VisSampleDataPrep(config).run();
    }
}
