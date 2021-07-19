package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.IM_FILE_ID;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleIdsFile;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.getGeneStatus;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadPredictionData;
import static com.hartwig.hmftools.neo.cohort.PredictionData.expandHomozygous;
import static com.hartwig.hmftools.neo.cohort.StatusResults.NORMAL;
import static com.hartwig.hmftools.neo.cohort.StatusResults.SIM_TUMOR;
import static com.hartwig.hmftools.neo.cohort.StatusResults.STATUS_MAX;
import static com.hartwig.hmftools.neo.cohort.StatusResults.TUMOR;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.neo.PeptideData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.lucene.search.spans.SpanMultiTermQueryWrapper;
import org.jetbrains.annotations.NotNull;

public class PeptidePredictionAnalyser
{
    private final String mOutputDir;
    private final String mNeoDataDir;
    private final String mPredictionsDataDir;
    private final List<String> mSampleIds;
    private final List<String> mSpecificAlleles;

    private BufferedWriter mPeptideWriter;

    private final int mPeptideAffinityThreshold;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String PREDICTION_DATA_DIR = "prediction_data_dir";

    private static final String PEPTIDE_AFF_THRESHOLD = "pep_aff_threshold";
    private static final String SPECIFIC_ALLELES = "specific_alleles";

    public PeptidePredictionAnalyser(final CommandLine cmd)
    {
        mSampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), mSampleIds);

        mNeoDataDir = cmd.getOptionValue(NEO_DATA_DIR);
        mPredictionsDataDir = cmd.getOptionValue(PREDICTION_DATA_DIR);

        mPeptideAffinityThreshold = Integer.parseInt(cmd.getOptionValue(PEPTIDE_AFF_THRESHOLD, "0"));
        mSpecificAlleles = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_ALLELES))
        {
            Arrays.stream(cmd.getOptionValue(SPECIFIC_ALLELES).split(ITEM_DELIM, -1)).forEach(x -> mSpecificAlleles.add(x));
            NE_LOGGER.info("filtering for {} alleles: {}", mSpecificAlleles.size(), mSpecificAlleles);
        }

        mOutputDir = parseOutputDir(cmd);
        mPeptideWriter = null;
        initialisePeptideWriter();
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
            return;

        NE_LOGGER.info("processing {} samples", mSampleIds.size());

        // check required inputs and config
        int processed = 0;

        for(final String sampleId : mSampleIds)
        {
            processSample(sampleId);
            ++processed;

            if(processed > 0 && (processed % 100) == 0)
            {
                NE_LOGGER.info("processed {} samples", processed);
            }
        }

        closeBufferedWriter(mPeptideWriter);
    }

    private void processSample(final String sampleId)
    {
        List<PredictionData> predictions = loadPredictionData(sampleId, mPredictionsDataDir);
        Map<Integer,NeoEpitopeData> neoDataMap = loadNeoEpitopes(sampleId);

        // organise into map by peptide, avoiding repeated peptides
        for(int i = 0; i < predictions.size(); ++i)
        {
            PredictionData predData = predictions.get(i);

            if(!mSpecificAlleles.isEmpty() && !mSpecificAlleles.contains(predData.Allele))
                continue;

            if(predData.Affinity > mPeptideAffinityThreshold)
                continue;

            NeoEpitopeData neoData = neoDataMap.get(predData.NeId);

            if(neoData == null)
            {
                NE_LOGGER.warn("sample({}) neId({}) neo-data not found", sampleId, predData.NeId);
                continue;
            }

            writePeptideData(sampleId, predData, neoData);
        }
    }

    private Map<Integer,NeoEpitopeData> loadNeoEpitopes(final String sampleId)
    {
        Map<Integer,NeoEpitopeData> neoDataMap = Maps.newHashMap();

        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(mNeoDataDir, sampleId);

            final List<String> lines = Files.readAllLines(new File(neoEpitopeFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int neIdIndex = fieldsIndexMap.get("NeId");
            int varTypeIndex = fieldsIndexMap.get("VariantType");
            int varInfoIndex = fieldsIndexMap.get("VariantInfo");
            int geneNameUpIndex = fieldsIndexMap.get("GeneNameUp");
            int geneNameDownIndex = fieldsIndexMap.get("GeneNameDown");

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                int neId = Integer.parseInt(items[neIdIndex]);

                String geneNameUp = items[geneNameUpIndex];
                String geneNameDown = items[geneNameDownIndex];
                String geneName = geneNameUp.equals(geneNameDown) ? geneNameUp : geneNameUp + "_" + geneNameDown;

                neoDataMap.put(neId, new NeoEpitopeData(
                        neId, NeoEpitopeType.valueOf(items[varTypeIndex]), items[varInfoIndex], geneName));
            }

            NE_LOGGER.debug("sample({}) loaded {} neo-epitopes", sampleId, lines.size());
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) neo-epitope file: {}", sampleId, exception.toString());
        }

        return neoDataMap;
    }


    private void initialisePeptideWriter()
    {
        try
        {
            final String outputFileName = String.format("%sNEO_PREDICTIONS_%d.csv", mOutputDir, mPeptideAffinityThreshold);

            mPeptideWriter = createBufferedWriter(outputFileName, false);
            mPeptideWriter.write("SampleId,NeId,VarType,VarInfo,Genes");
            mPeptideWriter.write(",Peptide,Allele,Affinity,Presentation");
            mPeptideWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create peptide writer: {}", e.toString());
        }
    }

    private void writePeptideData(final String sampleId, final PredictionData prediction, final NeoEpitopeData neoData)
    {
        try
        {
            mPeptideWriter.write(String.format("%s,%d,%s,%s,%s",
                    sampleId, neoData.Id, neoData.VariantType, neoData.VariantInfo, neoData.GeneName));

            mPeptideWriter.write(String.format(",%s,%s,%.1f,%.4f",
                    prediction.Peptide, prediction.Allele, prediction.Affinity, prediction.Presentation));

            mPeptideWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neoepitope prediction data: {}", e.toString());
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(PEPTIDE_AFF_THRESHOLD, true, "Only write peptides with affinity less than this if > 0");
        options.addOption(SPECIFIC_ALLELES, true, "Specific alleles to filter for, separated by ';'");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        PeptidePredictionAnalyser.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        PeptidePredictionAnalyser samplePeptidePredictions = new PeptidePredictionAnalyser(cmd);
        samplePeptidePredictions.run();

        NE_LOGGER.info("cohort peptide predictions complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
