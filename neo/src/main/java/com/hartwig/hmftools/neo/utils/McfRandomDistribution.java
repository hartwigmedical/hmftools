package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRED_AFFINITY;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PRES_SCORE;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.TrainConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.RandomDistributionTask.generateDistributionBuckets;
import static com.hartwig.hmftools.neo.bind.RandomPeptideDistribution.initialiseWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.VectorUtils;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.TrainConfig;
import com.hartwig.hmftools.neo.bind.RandomDistributionTask;
import com.hartwig.hmftools.neo.bind.RandomPeptideConfig;
import com.hartwig.hmftools.neo.bind.RandomPeptideDistribution;
import com.hartwig.hmftools.neo.bind.ScoreDistributionData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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

            String distributionFilename = TrainConfig.formFilename("mcf_random_peptide_dist", mConfig.OutputDir, mConfig.OutputId);
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

        final Map<String,Integer> otherColumns = Maps.newLinkedHashMap();

        if(!loadBindData(filename, Lists.newArrayList(), allelePeptideData, otherColumns))
            return;

        if(!mRandomDistribution.loadData())
            return;

        try
        {
            String peptideFilename = TrainConfig.formFilename("mcf_validation_peptide_scores", mConfig.OutputDir, mConfig.OutputId);

            BufferedWriter peptideWriter = createBufferedWriter(peptideFilename, false);
            peptideWriter.write("Allele,Peptide,RankPerc,PredictedAffinity,AffinityPerc,PresentationScore,PresentationPerc");
            peptideWriter.newLine();

            int predAffinityIndex = otherColumns.get(FLD_PRED_AFFINITY);
            int affinityPercIndex = otherColumns.get("AffinityPercentile");
            int presScoreIndex = otherColumns.get(FLD_PRES_SCORE);
            int presPercIndex = otherColumns.get("PresentationPercentile");

            for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : allelePeptideData.entrySet())
            {
                final String allele = alleleEntry.getKey();

                final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

                for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
                {
                    List<BindData> bindDataList = pepLenEntry.getValue();

                    for(BindData bindData : bindDataList)
                    {
                        double presentationScore = Double.parseDouble(bindData.getOtherData().get(presScoreIndex));
                        double presentationPercentile = Double.parseDouble(bindData.getOtherData().get(presPercIndex));
                        double predictedAffinity = Double.parseDouble(bindData.getOtherData().get(predAffinityIndex));
                        double affinityPerc = Double.parseDouble(bindData.getOtherData().get(affinityPercIndex));
                        double scoreValue = mUsePresentation ? presentationScore : predictedAffinity;
                        double scoreRank = mRandomDistribution.getScoreRank(allele, bindData.peptideLength(), scoreValue);

                        peptideWriter.write(String.format("%s,%s,%.6f,%.2f,%.6f,%.6f,%.6f",
                                allele, bindData.Peptide, scoreRank, predictedAffinity,
                                affinityPerc, presentationScore, presentationPercentile));
                        peptideWriter.newLine();
                    }
                }
            }

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

        // Allele,Peptide,PredictedAffinity,PresentationScore
        final Map<String,Map<Integer,List<BindData>>> allelePeptideData = Maps.newHashMap();

        final Map<String,Integer> otherColumns = Maps.newLinkedHashMap();

        if(!loadBindData(filename, Lists.newArrayList(), allelePeptideData, otherColumns))
            return;

        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : allelePeptideData.entrySet())
        {
            List<BindData> alleleBindList = Lists.newArrayList();
            alleleEntry.getValue().values().forEach(x -> alleleBindList.addAll(x));
            generateDistribution(alleleEntry.getKey(), alleleBindList, otherColumns);
        }

        return;
    }

    private void generateDistribution(final String allele, final List<BindData> bindDataList, final Map<String,Integer> otherColumns)
    {
        NE_LOGGER.info("allele({}) building distribution with {} random peptide predictions", allele, bindDataList.size());

        Map<Integer,List<Double>> pepLenScoresMap = Maps.newHashMap();

        int predAffinityIndex = otherColumns.get(FLD_PRED_AFFINITY);
        int presScoreIndex = otherColumns.get(FLD_PRES_SCORE);

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
            double presentationScore = Double.parseDouble(bindData.getOtherData().get(presScoreIndex));
            double predictedAffinity = Double.parseDouble(bindData.getOtherData().get(predAffinityIndex));

            double scoreValue = mUsePresentation ? presentationScore : predictedAffinity;
            boolean isAscending = mUsePresentation ? false : true;
            VectorUtils.optimisedAdd(peptideScores, scoreValue, isAscending);

            ++count;

            if(count > 0 && (count % 500000) == 0)
            {
                NE_LOGGER.debug("added {} sorted random peptide scores", count);
            }
        }

        final List<double[]> discreteScoreData = generateDistributionBuckets();

        for(Map.Entry<Integer,List<Double>> entry : pepLenScoresMap.entrySet())
        {
            NE_LOGGER.debug("allele({}) peptideLength({}) writing distribution", allele, entry.getKey());

            List<ScoreDistributionData> scoreDistribution = RandomDistributionTask.generateDistribution(
                    allele, entry.getKey(), entry.getValue(), discreteScoreData);

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
