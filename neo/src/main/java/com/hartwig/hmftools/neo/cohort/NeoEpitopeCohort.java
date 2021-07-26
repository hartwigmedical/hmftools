package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.VAR_INFO_DELIM;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.INFRAME_DELETION;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.neo.NeoEpitopeType.MISSENSE;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.IM_FILE_ID;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleIdsFile;
import static com.hartwig.hmftools.neo.predict.NeoPredictionsConfig.PREDICTIONS_FILE_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class NeoEpitopeCohort
{
    private final String mOutputDir;
    private final String mNeoDataDir;
    private final String mPredictionsDataDir;
    private final List<String> mSampleIds;

    private BufferedWriter mCohortWriter;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String PREDICTION_DATA_DIR = "prediction_data_dir";

    private static final double PREDICTION_FACTOR = 2;

    public NeoEpitopeCohort(final CommandLine cmd)
    {
        mSampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), mSampleIds);

        mNeoDataDir = cmd.getOptionValue(NEO_DATA_DIR);
        mPredictionsDataDir = cmd.getOptionValue(PREDICTION_DATA_DIR);
        mOutputDir = parseOutputDir(cmd);
        mCohortWriter = null;
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
            return;

        initialiseWriter();

        NE_LOGGER.info("processing {} samples", mSampleIds.size());

        // check required inputs and config
        int processed = 0;

        for(final String sampleId : mSampleIds)
        {
            processSampleNeoEpitopes(sampleId);
            ++processed;

            if(processed > 0 && (processed % 100) == 0)
            {
                NE_LOGGER.info("processed {} samples", processed);
            }
        }

        closeBufferedWriter(mCohortWriter);
    }

    private void processSampleNeoEpitopes(final String sampleId)
    {
        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(mNeoDataDir, sampleId);
            String rnaNeoEpitopeFile = neoEpitopeFile.replace(IM_FILE_ID, ISF_FILE_ID);

            boolean hasRnaData = Files.exists(Paths.get(rnaNeoEpitopeFile));

            final List<String> lines = hasRnaData ?
                    Files.readAllLines(new File(rnaNeoEpitopeFile).toPath()) : Files.readAllLines(new File(neoEpitopeFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            Integer fragIndex = fieldsIndexMap.get("FragmentsNovel");
            Integer baseDepthUp = fieldsIndexMap.get("BaseDepthUp");
            Integer baseDepthDown = fieldsIndexMap.get("BaseDepthDown");

            final Map<Integer,double[]> predictionSummary = loadPredictionSummary(sampleId);

            for(String line : lines)
            {
                NeoEpitopeFile neoData = null;

                try
                {
                    neoData = NeoEpitopeFile.fromString(line, false);
                }
                catch(Exception e)
                {
                    NE_LOGGER.error("sample({}) error({}) parsing neo data: {}", sampleId, e.toString(), line);
                    continue;
                }

                int rnaFragCount = 0;
                int[] rnaBaseDepth = {0, 0};

                if(hasRnaData)
                {
                    final String[] items = line.split(DELIMITER, -1);
                    rnaFragCount = Integer.parseInt(items[fragIndex]);
                    rnaBaseDepth[SE_START] =  Integer.parseInt(items[baseDepthUp]);
                    rnaBaseDepth[SE_END] =  Integer.parseInt(items[baseDepthDown]);
                }

                final double[] peptideValues = predictionSummary.get(neoData.Id);

                writeData(sampleId, neoData, hasRnaData, rnaFragCount, rnaBaseDepth, peptideValues);
            }

            NE_LOGGER.debug("sample({}) loaded {} neo-epitopes {}", sampleId, lines.size(), hasRnaData ? "with RNA" : "");
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) neo-epitope file: {}", sampleId, exception.toString());
        }
    }

    private static final int NINE_MER_COUNT = 0;
    private static final int MIN_AFFINITY = 1;
    private static final int MAX_PRESENTATION = 2;

    private Map<Integer,double[]> loadPredictionSummary(final String sampleId)
    {
        final Map<Integer,double[]> predictionSummary = Maps.newHashMap();

        try
        {
            final String predictionsFile = mPredictionsDataDir + sampleId + PREDICTIONS_FILE_ID;
            final List<String> fileContents = Files.readAllLines(new File(predictionsFile).toPath());

            if(fileContents.isEmpty())
                return predictionSummary;

            final String header = fileContents.get(0);
            fileContents.remove(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            // NeId,HlaAllele,Peptide,affinity,affinity_percentile,processing_score,presentation_score,presentation_percentile
            int neIdIndex = fieldsIndexMap.get("NeId");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int affinityIndex = fieldsIndexMap.get("affinity");
            int presentationIndex = fieldsIndexMap.get("presentation_score");

            int currentNeId = -1;
            double[] predictionValues = null;
            Set<String> ninemerPeptides = Sets.newHashSet();

            for(String data : fileContents)
            {
                final String[] items = data.split(DELIMITER);

                int neId = Integer.parseInt(items[neIdIndex]);
                double affinity = Double.parseDouble(items[affinityIndex]);
                double presentation = Double.parseDouble(items[presentationIndex]);
                String peptide = items[peptideIndex];

                if(currentNeId != neId)
                {
                    currentNeId = neId;
                    predictionValues = new double[MAX_PRESENTATION + 1];
                    predictionValues[MIN_AFFINITY] = affinity;
                    predictionSummary.put(neId, predictionValues);
                    ninemerPeptides.clear();
                }

                predictionValues[MAX_PRESENTATION] = max(predictionValues[MAX_PRESENTATION], presentation);

                if(affinity < predictionValues[MIN_AFFINITY])
                    predictionValues[MIN_AFFINITY] = affinity;

                if(peptide.length() == 9)
                {
                    ninemerPeptides.add(peptide);
                    predictionValues[NINE_MER_COUNT] = ninemerPeptides.size();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load sample predictions: {}", e.toString());
        }

        return predictionSummary;
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "IMU_NEO_EPITOPES.csv";

            mCohortWriter = createBufferedWriter(outputFileName, false);
            mCohortWriter.write("SampleId,NeId,VariantType,VariantInfo,JunctionCopyNumber");
            mCohortWriter.write(",GeneIdUp,GeneIdDown,GeneNameUp,GeneNameDown");
            mCohortWriter.write(",NmdMin,NmdMax,CodingBasesLengthMin,CodingBasesLengthMax,FusedIntronLength");
            mCohortWriter.write(",SkippedDonors,SkippedAcceptors,UpTranscripts,DownTranscripts");
            mCohortWriter.write(",TpmCancerUp,TpmCohortUp,TpmCancerDown,TpmCohortDown");
            mCohortWriter.write(",HasRna,RnaFragCount,RnaDepthUp,RnaDepthDown");
            mCohortWriter.write(",NineMers,MinAffinity,MaxPresentation");
            mCohortWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create neo-epitope writer: {}", e.toString());
        }
    }

    private void writeData(
            final String sampleId, final NeoEpitopeFile neData,
            boolean hasRna, int rnaFragCount, final int[] rnaBaseDepth, final double[] peptideValues)
    {
        try
        {
            // correction for missense which are inframe INDELs
            NeoEpitopeType variantType = neData.VariantType;

            if(variantType == MISSENSE)
            {
                String[] variantInfo = neData.VariantInfo.split(VAR_INFO_DELIM, -1);
                if(variantInfo.length == 4 && variantInfo[2].length() != variantInfo[3].length())
                {
                    variantType = variantInfo[3].length() > variantInfo[2].length() ? INFRAME_INSERTION : INFRAME_DELETION;
                }
            }

            mCohortWriter.write(String.format("%s,%d,%s,%s,%.1f",
                    sampleId, neData.Id, variantType, neData.VariantInfo, neData.CopyNumber));

            mCohortWriter.write(String.format(",%s,%s,%s,%s",
                    neData.GeneIds[FS_UP], neData.GeneIds[FS_DOWN], neData.GeneNames[FS_UP], neData.GeneNames[FS_DOWN]));

            mCohortWriter.write(String.format(",%d,%d,%d,%d,%d",
                    neData.NmdBases[0], neData.NmdBases[1], neData.CodingBasesLength[0], neData.CodingBasesLength[1], neData.FusedIntronLength));

            mCohortWriter.write(String.format(",%d,%d,%s,%s",
                    neData.SkippedAcceptorsDonors[FS_UP], neData.SkippedAcceptorsDonors[FS_DOWN],
                    neData.Transcripts[FS_UP], neData.Transcripts[FS_DOWN]));

            mCohortWriter.write(String.format(",%6.3e,%6.3e,%6.3e,%6.3e",
                    neData.CancerTpmTotal[FS_UP], neData.CohortTpmTotal[FS_UP], neData.CancerTpmTotal[FS_DOWN], neData.CohortTpmTotal[FS_DOWN]));

            mCohortWriter.write(String.format(",%s,%d,%d,%d",
                    hasRna, rnaFragCount, rnaBaseDepth[SE_START], rnaBaseDepth[SE_END]));

            if(peptideValues != null)
            {
                mCohortWriter.write(String.format(",%.0f,%.2f,%.6f",
                        peptideValues[NINE_MER_COUNT], peptideValues[MIN_AFFINITY], peptideValues[MAX_PRESENTATION]));
            }
            else
            {
                mCohortWriter.write(",0,0,0");
            }

            mCohortWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        NeoEpitopeCohort.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NeoEpitopeCohort neoEpitopeCohort = new NeoEpitopeCohort(cmd);
        neoEpitopeCohort.run();

        NE_LOGGER.info("Neo-epitope annotations complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
