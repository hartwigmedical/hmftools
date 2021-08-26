package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRED_AFFINITY;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRES_SCORE;
import static com.hartwig.hmftools.neo.bind.BindCommon.RANDOM_SOURCE;
import static com.hartwig.hmftools.neo.bind.BindData.cleanAllele;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.RandomPeptideDistribution.initialiseWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.VectorUtils;
import com.hartwig.hmftools.common.stats.AucCalc;
import com.hartwig.hmftools.common.stats.AucData;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BinderConfig;
import com.hartwig.hmftools.neo.bind.RandomPeptideConfig;
import com.hartwig.hmftools.neo.bind.RandomPeptideDistribution;
import com.hartwig.hmftools.neo.bind.ScoreDistributionData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

// rather than use/trust McfFlurry's affinity and presentation score percentiles, rebuild a distribution for each of them using their
// scores on a set of random peptides (eg 100K)
public class McfRandomDistribution
{
    private final RandomPeptideConfig mConfig;
    private final String mMcfPredictionsFile;
    private final String mValidationDataFile;
    private final boolean mUsePresentation;

    private final RandomPeptideDistribution mRandomDistribution;
    private BufferedWriter mDistributionWriter;

    private static final String PREDICTIONS_FILE = "mcf_rand_predictions_file";
    private static final String VALIDATION_FILE = "validation_data_file";
    private static final String USE_PRESENTATION = "use_presentation";

    public McfRandomDistribution(final CommandLine cmd)
    {
        mMcfPredictionsFile = cmd.getOptionValue(PREDICTIONS_FILE);
        mValidationDataFile = cmd.getOptionValue(VALIDATION_FILE);

        mConfig = new RandomPeptideConfig(cmd);
        mUsePresentation = cmd.hasOption(USE_PRESENTATION);

        mRandomDistribution = new RandomPeptideDistribution(mConfig);

        mDistributionWriter = null;
    }

    public void run()
    {
        if(mMcfPredictionsFile != null)
        {
            NE_LOGGER.info("loading MCF random peptide predictions from file({}) using {}",
                    mMcfPredictionsFile, mUsePresentation ? "presentation" : "affinity");

            String distributionFilename = BinderConfig.formFilename("mcf_random_peptide_dist", mConfig.OutputDir, mConfig.OutputId);
            mDistributionWriter = initialiseWriter(distributionFilename);

            processFile(mMcfPredictionsFile);
            closeBufferedWriter(mDistributionWriter);

            NE_LOGGER.info("MCF random peptide affinity distribution build complete");
            return;
        }

        if(mValidationDataFile != null)
        {
            processValidationData(mValidationDataFile);
        }
    }

    private void processValidationData(final String filename)
    {
        final Map<String,Map<Integer,List<BindData>>> allelePeptideData = Maps.newHashMap();

        if(!loadBindData(filename, true, Lists.newArrayList(), Lists.newArrayList(), allelePeptideData))
            return;

        if(!mRandomDistribution.loadData())
            return;

        try
        {
            String alleleFilename = BinderConfig.formFilename("mcf_validation_summary", mConfig.OutputDir, mConfig.OutputId);

            NE_LOGGER.info("writing score results to {}", alleleFilename);

            BufferedWriter alleleWriter = createBufferedWriter(alleleFilename, false);
            alleleWriter.write("Allele,PeptideLength,DataCount,AUC");
            alleleWriter.newLine();

            String peptideFilename = BinderConfig.formFilename("mcf_validation_peptide_scores", mConfig.OutputDir, mConfig.OutputId);

            BufferedWriter peptideWriter = createBufferedWriter(peptideFilename, false);
            peptideWriter.write("Allele,Peptide,RankPerc,PredictedAffinity,AffinityPerc,PresentationScore,PresentationPerc");
            peptideWriter.newLine();

            for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : allelePeptideData.entrySet())
            {
                final String allele = alleleEntry.getKey();

                final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

                for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
                {
                    int peptideLength = pepLenEntry.getKey();
                    List<BindData> bindDataList = pepLenEntry.getValue();

                    List<AucData> aucData = Lists.newArrayList();

                    for(BindData bindData : bindDataList)
                    {
                        double scoreValue = mUsePresentation ? bindData.presentationScore() : bindData.predictedAffinity();
                        double scoreRank = mRandomDistribution.getScoreRank(allele, bindData.peptideLength(), scoreValue);

                        peptideWriter.write(String.format("%s,%s,%.6f,%.2f,%.6f,%.6f,%.6f",
                                allele, bindData.Peptide, scoreRank, bindData.predictedAffinity(),
                                bindData.affinityPercentile(), bindData.presentationScore(), bindData.presentationPercentile()));
                        peptideWriter.newLine();

                        aucData.add(new AucData(true, scoreRank, true));
                    }

                    double aucPerc = AucCalc.calcPercentilesAuc(aucData, Level.TRACE);

                    NE_LOGGER.info(String.format("allele(%s) peptideLength(%d) items(%d) AUC(%.4f)",
                            allele, peptideLength, bindDataList.size(), aucPerc));

                    alleleWriter.write(String.format("%s,%d,%d,%.4f", allele, peptideLength, bindDataList.size(), aucPerc));
                    alleleWriter.newLine();
                }
            }

            alleleWriter.close();
            peptideWriter.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to init allele summary writer: {}", e.toString());
            return;
        }
    }

    private void processFile(final String filename)
    {
        // read the predictions for each allele in turn and create a distribution from it's affinity scores
        if(filename == null || !Files.exists(Paths.get(filename)))
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int affinityIndex = fieldsIndexMap.get(FLD_PRED_AFFINITY);
            Integer presentationIndex = fieldsIndexMap.get(FLD_PRES_SCORE);

            String currentAllele = "";
            List<BindData> bindDataList = null;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = cleanAllele(items[alleleIndex]);
                String peptide = items[peptideIndex];
                double predictedAffinity = Double.parseDouble(items[affinityIndex]);
                double presScore = presentationIndex != null ? Double.parseDouble(items[presentationIndex]) : -1;

                BindData bindData = new BindData(allele, peptide, RANDOM_SOURCE);
                bindData.setPredictionData(predictedAffinity, -1, presScore, -1);

                if(!currentAllele.equals(allele))
                {
                    if(bindDataList != null)
                        generateDistribution(currentAllele, bindDataList);

                    currentAllele = allele;
                    bindDataList = Lists.newArrayList();
                }

                bindDataList.add(bindData);
            }

            generateDistribution(currentAllele, bindDataList);

        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read MCF random predictions data file: {}", e.toString());
            return;
        }

        return;
    }

    private void generateDistribution(final String allele, final List<BindData> bindDataList)
    {
        NE_LOGGER.info("allele({}) building distribution with {} random peptide predictions", allele, bindDataList.size());

        Map<Integer,List<Double>> pepLenScoresMap = Maps.newHashMap();

        int count = 0;

        for(BindData bindData : bindDataList)
        {
            List<Double> peptideScores = pepLenScoresMap.get(bindData.peptideLength());

            if(peptideScores == null)
            {
                peptideScores = Lists.newArrayList();
                pepLenScoresMap.put(bindData.peptideLength(), peptideScores);
            }

            // affinity goes from low to high as binding gets worse, presentation is the opposite
            double scoreValue = mUsePresentation ? bindData.presentationScore() : bindData.predictedAffinity();
            boolean isAscending = mUsePresentation ? false : true;
            VectorUtils.optimisedAdd(peptideScores, scoreValue, isAscending);

            ++count;

            if(count > 0 && (count % 500000) == 0)
            {
                NE_LOGGER.debug("added {} sorted random peptide scores", count);
            }
        }

        for(Map.Entry<Integer,List<Double>> entry : pepLenScoresMap.entrySet())
        {
            NE_LOGGER.debug("allele({}) peptideLength({}) writing distribution", allele, entry.getKey());

            List<ScoreDistributionData> scoreDistribution = mRandomDistribution.generateDistribution(
                    allele, entry.getKey(), entry.getValue());

            RandomPeptideDistribution.writeDistribution(mDistributionWriter, scoreDistribution);
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(PREDICTIONS_FILE, true, "MCF predictions file");
        options.addOption(VALIDATION_FILE, true, "Binding validation file");
        options.addOption(USE_PRESENTATION, false, "Rank and score using presentation instead of affinity");
        addLoggingOptions(options);
        addOutputDir(options);
        options.addOption(OUTPUT_ID, true, "Output file id");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        McfRandomDistribution mcfRandomDistribution = new McfRandomDistribution(cmd);
        mcfRandomDistribution.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
