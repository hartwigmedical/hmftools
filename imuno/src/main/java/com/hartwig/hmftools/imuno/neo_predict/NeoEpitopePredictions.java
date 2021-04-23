package com.hartwig.hmftools.imuno.neo_predict;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.imuno.neo.NeoConfig.NEO_EPITOPE_FILE_ID;
import static com.hartwig.hmftools.imuno.neo_predict.NeoPredictionsConfig.PREDICTIONS_FILE_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.imuno.neo.SampleData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class NeoEpitopePredictions
{
    private final NeoPredictionsConfig mConfig;

    private SampleData mCurrentSample;
    private BufferedWriter mWriter;

    public NeoEpitopePredictions(final CommandLine cmd)
    {
        mConfig = new NeoPredictionsConfig(cmd);
        mWriter = null;
        mCurrentSample = null;
    }

    public void run()
    {
        if(mConfig.Samples.isEmpty())
            return;

        if(mConfig.Samples.size() == 1)
        {
            IM_LOGGER.info("processing sample({})", mConfig.Samples.get(0).Id);
        }
        else
        {
            IM_LOGGER.info("processing {} samples", mConfig.Samples.size());
        }

        // check required inputs and config
        int sampleCount = 0;
        int nextLog = 100;
        for(final SampleData sample : mConfig.Samples)
        {
            mCurrentSample = sample;


            final List<NeoPredictionData> predictions = loadPredictionData(sample.Id);
            final List<NeoEpitopeFile> neoEpitopes = loadNeoEpitopeData(sample.Id);

            IM_LOGGER.debug("sample({}) loaded {} fusions and {} point mutations",
                    sample.Id, predictions.size());

            if(!mConfig.WriteCohortFile)
            {
                closeBufferedWriter(mWriter);
                mWriter = null;
            }
        }

        if(mConfig.WriteCohortFile)
        {
            closeBufferedWriter(mWriter);
        }
    }

    private final List<NeoPredictionData> loadPredictionData(final String sampleId)
    {
        final List<NeoPredictionData> predictions = Lists.newArrayList();

        try
        {
            final String samplePredictionsFile = mConfig.SampleDataDir + sampleId + PREDICTIONS_FILE_ID;
            final List<String> fileContents = Files.readAllLines(new File(samplePredictionsFile).toPath());

            if(fileContents.isEmpty())
                return predictions;

            final String header = fileContents.get(0);
            fileContents.remove(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            // NeId,HlaAllele,Peptide,affinity,affinity_percentile,processing_score,presentation_score,presentation_percentile
            int neIdIndex = fieldsIndexMap.get("NeId");
            int hlaIndex = fieldsIndexMap.get("HlaAllele");
            int peptideIndex = fieldsIndexMap.get("Peptide");

            NeoPredictionData nedPredData = null;

            for(String data : fileContents)
            {
                final String[] items = data.split(DELIMITER);

                int neId = Integer.parseInt(items[neIdIndex]);

                if(nedPredData == null || nedPredData.NeId != neId)
                {
                    nedPredData = new NeoPredictionData(neId);
                    predictions.add(nedPredData);
                }

                nedPredData.Predictions.add(new HlaPeptidePrediction(
                        items[hlaIndex], items[peptideIndex], PredictionValues.fromCsv(items, peptideIndex+1)));
            }
        }
        catch (IOException e)
        {
            IM_LOGGER.warn("failed to load sample predictions: {}", e.toString());
        }

        return predictions;
    }

    private void writeData()
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                String outputFileName = mConfig.OutputDir;

                if(mConfig.WriteCohortFile)
                    outputFileName += "IMU_NEO_EPITOPES.csv";
                else
                    outputFileName += mCurrentSample.Id + ".imu.neo_epitopes.csv";

                mWriter = createBufferedWriter(outputFileName, false);

                mWriter.newLine();
            }

            if(mConfig.WriteCohortFile)
                mWriter.write(String.format("%s,", mCurrentSample.Id));

            // final NeoEpitopeFile neFile = neData.toFile(neId, upTransNames, downTransNames, tpmCancer, tpmCohort);
            // mWriter.write(NeoEpitopeFile.toString(neFile));
            mWriter.newLine();
        }
        catch (final IOException e)
        {
            IM_LOGGER.error("error writing neo-epitope output file: {}", e.toString());
        }
    }

    private List<NeoEpitopeFile> loadNeoEpitopeData(final String sampleId)
    {
        final List<NeoEpitopeFile> neDataList = Lists.newArrayList();

        try
        {
            final String sampleNeDataFile = mConfig.SampleDataDir + sampleId + NEO_EPITOPE_FILE_ID;
            neDataList.addAll(NeoEpitopeFile.read(sampleNeDataFile));
        }
        catch(IOException exception)
        {
            IM_LOGGER.error("failed to read neo-epitope file:", exception.toString());
        }

        return neDataList;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoPredictionsConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        NeoEpitopePredictions neoEpitopePredictions = new NeoEpitopePredictions(cmd);
        neoEpitopePredictions.run();

        IM_LOGGER.info("Neo-epitope prediction analysis complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
