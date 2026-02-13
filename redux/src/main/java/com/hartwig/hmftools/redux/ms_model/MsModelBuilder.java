package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.COEFFICIENTS;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.ERROR_RATES;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.VALIDATION;
import static com.hartwig.hmftools.redux.ms_model.MsModelParams.DEFAULT_MODEL_PARAMS;

import java.io.IOException;
import java.util.Collection;
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

        mModelParams = DEFAULT_MODEL_PARAMS;

        mSamplePurities = Maps.newHashMap();
        mSampleJitterCounts = Maps.newHashMap();
    }

    public void run()
    {
        RD_LOGGER.info("running MSI model builder for {} samples", mConfig.SampleIds.size());

        loadSampleFiles();

        MsModelCalculator modelCalculator = new MsModelCalculator(mModelParams);

        if(mConfig.Routines.contains(ERROR_RATES))
        {
            modelCalculator.calculateErrorRates(mSampleJitterCounts);

            // persist to file
            String modelCoefficientsFile = mConfig.ModelCoefficientsFile != null ?
                    mConfig.ModelCoefficientsFile : MsModelCoefficients.generateFilename(mConfig.OutputDir, mConfig.OutputId);
            MsModelCoefficients.write(modelCoefficientsFile, modelCalculator.getModelCoeffcients());
        }
        else
        {
            // load error rates
        }

        if(mConfig.Routines.contains(COEFFICIENTS))
        {
            modelCalculator.calculateCoefficients(mSamplePurities);
        }
        else
        {
            // load coefficients
        }

        if(mConfig.Routines.contains(VALIDATION))
        {
            modelCalculator.runValidation(mSamplePurities, mSampleJitterCounts);
        }

        RD_LOGGER.info("MSI model builder complete");
    }

    private void loadSampleFiles()
    {
        try
        {
            for(String sampleId : mConfig.SampleIds)
            {
                String purplePurityFile = PurplePurity.generateFilename(mConfig.PurpleDir, sampleId);
                PurplePurity purplePurity = PurplePurity.read(purplePurityFile);
                mSamplePurities.put(sampleId, purplePurity);

                String msSitesFile = JitterCountsTableFile.generateFilename(mConfig.ReduxDir, sampleId);

                Collection<JitterCountsTable> jitterCounts = JitterCountsTableFile.read(msSitesFile);
                mSampleJitterCounts.put(sampleId, jitterCounts);
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
