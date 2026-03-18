package com.hartwig.hmftools.redux.ms_model;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.COEFFICIENTS;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.ERROR_RATES;
import static com.hartwig.hmftools.redux.ms_model.MsModelConfig.TrainingRoutines.VALIDATION;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterCountsTableFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.Level;

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
            modelCalculator.loadErrorRates(mConfig.ModelErroRatesFile);
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
            modelCalculator.loadCoefficients(mConfig.ModelCoefficientsFile);
        }

        if(mConfig.Routines.contains(VALIDATION))
        {
            runValidation(modelCalculator, mSamplePurities, mSampleJitterCounts, mConfig);
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

    public void runValidation(
            final MsModelCalculator modelCalculator, final Map<String,PurplePurity> samplePurities,
            final Map<String,Collection<JitterCountsTable>> sampleJitterCounts, final MsModelConfig config)
    {
        RD_LOGGER.info("evaluating samples");

        Level logLevel = mSamplePurities.size() > 100 ? Level.TRACE : Level.DEBUG;

        try
        {
            String filename = config.OutputDir + File.separator + "ms_model_evaluation";
            String calcsFilename = config.OutputDir + File.separator + "ms_model_calcs";

            if(config.OutputId != null)
            {
                filename += "." + config.OutputId;
                calcsFilename += "." + config.OutputId;
            }

            filename += TSV_EXTENSION;
            calcsFilename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_SAMPLE_ID).add("Purity").add("MsIndelsPerMb").add("PredictedValue");
            writer.write(sj.toString());
            writer.newLine();

            BufferedWriter calcsWriter = createBufferedWriter(calcsFilename, false);
            sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_SAMPLE_ID).add("RepeatUnit").add("RepeatCount").add("AdjustedErrorRate").add("PredictedValue");
            calcsWriter.write(sj.toString());
            calcsWriter.newLine();

            for(Map.Entry<String,PurplePurity> entry : samplePurities.entrySet())
            {
                String sampleId = entry.getKey();

                PurplePurity purplePurity = entry.getValue();

                Collection<JitterCountsTable> jitterCounts = sampleJitterCounts.get(sampleId);

                List<RepeatUnitData> repeatUnitDataList = modelCalculator.buildRepeatUnitData(jitterCounts);

                List<Double> repeatUnitPredictedValues = modelCalculator.calcMsIndelPerMbValues(repeatUnitDataList);

                for(int i = 0; i < repeatUnitDataList.size(); ++i)
                {
                    RepeatUnitData repeatUnitData = repeatUnitDataList.get(i);
                    double predictedValue = repeatUnitPredictedValues.get(i);

                    sj = new StringJoiner(TSV_DELIM);
                    sj.add(sampleId);
                    sj.add(repeatUnitData.RepeatUnit);
                    sj.add(String.valueOf(repeatUnitData.repeatCount()));
                    sj.add(format("%4.3e", repeatUnitData.adjustedErrorRate()));
                    sj.add(format("%.4f", predictedValue));
                    calcsWriter.write(sj.toString());
                    calcsWriter.newLine();
                }

                double predictedMsIndelsPerMb = modelCalculator.calcMsIndelPerMb(repeatUnitPredictedValues);

                RD_LOGGER.log(logLevel, format("sample(%s) purity(%.3f) msIndelsPerMb(actual=%.4f predicted=%.4f)",
                        sampleId, purplePurity.Purity, purplePurity.MsIndelsPerMb, predictedMsIndelsPerMb));

                sj = new StringJoiner(TSV_DELIM);
                sj.add(sampleId);
                sj.add(format("%.4f", purplePurity.Purity));
                sj.add(format("%.4f", purplePurity.MsIndelsPerMb));
                sj.add(format("%.4f", predictedMsIndelsPerMb));
                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
            calcsWriter.close();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to create sample evaluation file: {}", e.toString());
            System.exit(1);
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
