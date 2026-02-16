package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.COEFFICIENTS;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.ERROR_RATES;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.VALIDATION;
import static com.hartwig.hmftools.redux.ms_model.MsModelParams.DEFAULT_MODEL_PARAMS;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterCountsTableFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MsModelBuilder
{
    private final MsModelConfig mConfig;
    private final Map<String,PurplePurity> mSamplePurities;
    private final Map<String,Collection<JitterCountsTable>> mSampleJitterCounts;

    private final MsModelParams mModelParams;

    public MsModelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new MsModelConfig(configBuilder);
        mModelParams = MsModelParams.fromConfig(configBuilder);

        mSamplePurities = Maps.newHashMap();
        mSampleJitterCounts = Maps.newHashMap();
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        RD_LOGGER.info("running MSI model builder for {} samples", mConfig.SampleIds.size());

        loadSampleFiles();

        MsModelCalculator modelCalculator = new MsModelCalculator(mModelParams);

        if(mConfig.Routines.contains(ERROR_RATES))
        {
            modelCalculator.calculateErrorRates(mSampleJitterCounts);

            // persist to file
            String modelErrorRatesFile = mConfig.ModelErroRatesFile != null ?
                    mConfig.ModelErroRatesFile : MsModelErrorRates.generateFilename(mConfig.OutputDir, mConfig.OutputId);
            MsModelErrorRates.write(modelErrorRatesFile, modelCalculator.getModelErrorRates());
        }
        else
        {
            // load error rates
            List<MsModelErrorRates> modelErrorRates = MsModelErrorRates.read(mConfig.ModelErroRatesFile);
            modelCalculator.loadModelErrorRates(modelErrorRates);
        }

        if(mConfig.Routines.contains(COEFFICIENTS))
        {
            modelCalculator.calculateCoefficients(mSamplePurities);

            String modelCoefficientsFile = mConfig.ModelCoefficientsFile != null ?
                    mConfig.ModelCoefficientsFile : MsModelCoefficients.generateFilename(mConfig.OutputDir, mConfig.OutputId);
            MsModelCoefficients.write(modelCoefficientsFile, modelCalculator.getModelCoeffcients());
        }
        else
        {
            // load coefficients
            List<MsModelCoefficients> modelCoefficients = MsModelCoefficients.read(mConfig.ModelCoefficientsFile);
            modelCalculator.loadModelCoeffcients(modelCoefficients);
        }

        if(mConfig.Routines.contains(VALIDATION))
        {
            modelCalculator.runValidation(mSamplePurities, mSampleJitterCounts, mConfig);
        }

        RD_LOGGER.info("MSI model builder, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void loadSampleFiles()
    {
        try
        {
            int count = 0;

            for(String sampleId : mConfig.SampleIds)
            {
                String purplePurityFile = PurplePurity.generateFilename(mConfig.PurpleDir, sampleId);
                String msSitesFile = JitterCountsTableFile.generateFilename(mConfig.ReduxDir, sampleId);

                if(!Files.exists(Paths.get(msSitesFile)))
                {
                    msSitesFile = msSitesFile.replace(".redux", "");
                }

                if(!Files.exists(Paths.get(msSitesFile)) || !Files.exists(Paths.get(purplePurityFile)))
                {
                    RD_LOGGER.warn("sample({}) missing files: present purity({}) msTable({})",
                            sampleId,
                            Files.exists(Paths.get(purplePurityFile)),
                            Files.exists(Paths.get(msSitesFile)));
                    continue;
                }

                PurplePurity purplePurity = PurplePurity.read(purplePurityFile);
                mSamplePurities.put(sampleId, purplePurity);

                Collection<JitterCountsTable> jitterCounts = JitterCountsTableFile.read(msSitesFile);
                mSampleJitterCounts.put(sampleId, jitterCounts);

                ++count;

                if((count % 1000) == 0)
                {
                    RD_LOGGER.debug("loaded {} sample files", count);
                }
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to load sample file: {}", e.toString());
        }
    }

    public static void main(final String... args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        MsModelConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MsModelBuilder msModelBuilder = new MsModelBuilder(configBuilder);
        msModelBuilder.run();
    }
}
