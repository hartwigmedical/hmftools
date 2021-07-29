package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindData.RANDOM_SOURCE;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.LIKELY_INCORRECT;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.TRAINING;
import static com.hartwig.hmftools.neo.bind.ScoringData.initMatrixWriter;
import static com.hartwig.hmftools.neo.bind.ScoringData.initPairDataWriter;
import static com.hartwig.hmftools.neo.bind.ScoringData.initFrequencyWriter;
import static com.hartwig.hmftools.neo.utils.AminoAcidFrequency.AMINO_ACID_FREQ_FILE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.stats.AucCalc;
import com.hartwig.hmftools.common.stats.AucData;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

public class BindTrainer
{
    private final BinderConfig mConfig;

    private final Map<String,List<BindData>> mAllelePeptideData;
    private final AminoAcidFrequency mAminoAcidFrequency;

    private final List<BindScoreMatrix> mBindMatrixList;
    private final Set<Integer> mDistinctPeptideLengths;

    public BindTrainer(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mAllelePeptideData = Maps.newHashMap();
        mDistinctPeptideLengths = Sets.newHashSet();

        mAminoAcidFrequency = new AminoAcidFrequency(cmd);
        mAminoAcidFrequency.loadFrequencies();

        mBindMatrixList = Lists.newArrayList();
    }

    public void run()
    {
        NE_LOGGER.info("running NeoBinder on {} alleles", mConfig.SpecificAlleles.isEmpty() ? "all" : mConfig.SpecificAlleles.size());

        if(!loadTrainingData() || !loadRandomPredictionsData())
        {
            System.exit(1);
        }

        if(mConfig.WriteScoreMatrix && mAminoAcidFrequency.getAminoAcidFrequencies().isEmpty())
        {
            NE_LOGGER.warn("no amino acid frequencies loaded");
        }

        processingBindingData();

        if(mConfig.RunScoring)
        {
            runScoring();
        }

        NE_LOGGER.info("NeoBinder complete");
    }

    private void processingBindingData()
    {
        if(mConfig.BindMatrixFile != null)
        {
            mBindMatrixList.addAll(BindScoreMatrix.loadFromCsv(mConfig.BindMatrixFile));
        }
        else
        {
            buildBindingScoreData();
        }
    }

    private void buildBindingScoreData()
    {
        BufferedWriter matrixWriter = mConfig.WriteScoreMatrix ?
                initMatrixWriter(mConfig.formFilename("score_matrix"), getMaxPeptideLength()) : null;

        BufferedWriter freqWriter = mConfig.WriteFrequencyData ?
                initFrequencyWriter(mConfig.formFilename("pos_frequency")) : null;

        BufferedWriter pairWriter = mConfig.CalcPairs ? initPairDataWriter(mConfig.formFilename("pair_score_prob")) : null;

        for(Map.Entry<String,List<BindData>> entry : mAllelePeptideData.entrySet())
        {
            final String allele = entry.getKey();
            List<BindData> peptideBindData = entry.getValue();

            NE_LOGGER.debug("allele({}) processing {} data items", allele, peptideBindData.size());

            Map<Integer,ScoringData> matrixMap = Maps.newHashMap();
            int currentLength = -1;
            ScoringData currentData = null;

            for(BindData bindData : peptideBindData)
            {
                int peptideLength = bindData.peptideLength();

                if(currentLength != peptideLength)
                {
                    currentLength = peptideLength;
                    currentData = matrixMap.get(peptideLength);

                    if(currentData == null)
                    {
                        currentData = new ScoringData(allele, peptideLength, mConfig.Constants);
                        matrixMap.put(peptideLength, currentData);
                    }
                }

                currentData.processBindData(bindData, mConfig.CalcPairs);
            }

            for(ScoringData scoringData : matrixMap.values())
            {
                // write results
                scoringData.logStats();

                BindScoreMatrix matrix = scoringData.createMatrix(mAminoAcidFrequency);
                mBindMatrixList.add(matrix);

                if(mConfig.WriteScoreMatrix)
                    scoringData.writeMatrixData(matrixWriter, matrix, getMaxPeptideLength());

                if(mConfig.WriteFrequencyData)
                    scoringData.writeFrequencyData(freqWriter);

                if(mConfig.CalcPairs)
                    scoringData.writePairData(allele, pairWriter);
            }
        }

        closeBufferedWriter(matrixWriter);
        closeBufferedWriter(freqWriter);
        closeBufferedWriter(pairWriter);
    }

    private void runScoring()
    {
        RandomPeptideDistribution randomDistribution = new RandomPeptideDistribution(mConfig);

        if(!randomDistribution.hasData())
        {
            randomDistribution.buildDistribution(mBindMatrixList, mAllelePeptideData);
        }

        NE_LOGGER.info("writing allele summaries");

        BufferedWriter alleleWriter = initAlleleSummaryWriter();

        // rank both the training and random peptide data using the newly created data and the random peptide distributions
        for(Map.Entry<String,List<BindData>> entry : mAllelePeptideData.entrySet())
        {
            String allele = entry.getKey();
            List<BindData> peptideBindData = entry.getValue();

            for(int peptideLength : mDistinctPeptideLengths)
            {
                BindScoreMatrix matrix = mBindMatrixList.stream()
                        .filter(x -> x.Allele.equals(allele) && x.PeptideCount == peptideLength).findFirst().orElse(null);

                for(BindData bindData : peptideBindData)
                {
                    if(bindData.peptideLength() != peptideLength)
                        continue;

                    bindData.Score = matrix.calcScore(bindData.Peptide);
                    bindData.RankPerc = randomDistribution.getScoreRank(allele, bindData.Score);
                }
            }

            writeAlleleSummary(alleleWriter, allele);
        }

        closeBufferedWriter(alleleWriter);

        NE_LOGGER.info("writing peptide scores");
        writePeptideScores();
    }

    private BufferedWriter initAlleleSummaryWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formFilename("allele_summary"), false);
            writer.write("Allele,TrainingCount,RandomCount,Auc,AucMcf");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to init allele summary writer: {}", e.toString());
            return null;
        }
    }

    private void writeAlleleSummary(final BufferedWriter writer, final String allele)
    {
        if(writer == null)
            return;

        try
        {
            final List<BindData> peptideBindData = mAllelePeptideData.get(allele);

            int trainingCount = 0;
            int randomCount = 0;

            List<AucData> alleleAucData = Lists.newArrayList();
            List<AucData> alleleAucMcfData = Lists.newArrayList();

            for(BindData bindData : peptideBindData)
            {

                if(bindData.isTraining())
                    ++trainingCount;
                else
                    ++randomCount;

                boolean isPositive = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;
                alleleAucData.add(new AucData(isPositive, bindData.Score));
                alleleAucMcfData.add(new AucData(isPositive, bindData.PredictedAffinity));
            }

            double auc = AucCalc.calcAuc(alleleAucData, Level.TRACE);
            double aucMcf = AucCalc.calcAuc(alleleAucMcfData, Level.TRACE);

            NE_LOGGER.info(String.format("allele(%s) peptides(train=%d, rand=%d) AUC(%.4f) mcfAUC(%.4f)",
                    allele, trainingCount, randomCount, auc, aucMcf));

            writer.write(String.format("%s,%d,%d,%.4f,%.4f", allele, trainingCount, randomCount, auc, aucMcf));
            writer.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write allele summary data: {}", e.toString());
        }
    }

    private void writePeptideScores()
    {
        if(mConfig.WritePeptideType == PeptideWriteType.NONE)
            return;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formFilename("peptide_scores"), false);
            writer.write("Allele,Peptide,Source,Score,Rank,Affinity,PredictedAffinity");
            writer.newLine();

            for(Map.Entry<String, List<BindData>> entry : mAllelePeptideData.entrySet())
            {
                String allele = entry.getKey();
                List<BindData> peptideBindData = entry.getValue();

                for(BindData bindData : peptideBindData)
                {
                    if(mConfig.WritePeptideType == LIKELY_INCORRECT)
                    {
                        boolean expectBinds = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;
                        boolean inTopPerc = bindData.RankPerc < 0.02;
                        if(expectBinds == inTopPerc)
                            continue;
                    }
                    else if(mConfig.WritePeptideType == TRAINING && !bindData.isTraining())
                    {
                        continue;
                    }

                    writer.write(String.format("%s,%s,%s,%.4f,%.4f,%.1f,%.1f",
                            allele, bindData.Peptide, bindData.Source,
                            bindData.Score, bindData.RankPerc, bindData.Affinity, bindData.PredictedAffinity));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide scores file: {}", e.toString());
        }
    }

    private boolean loadTrainingData()
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(mConfig.TrainingDataFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get("Allele");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int affinityIndex = fieldsIndexMap.get("Affinity");
            int predictedIndex = fieldsIndexMap.get("PredictedAffinity");
            int sourceIndex = fieldsIndexMap.get("Source");

            String currentAllele = "";
            List<BindData> currentBindList = null;

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = items[alleleIndex];

                if(!mConfig.SpecificAlleles.isEmpty() && !mConfig.SpecificAlleles.contains(allele))
                {
                    if(mConfig.SpecificAlleles.size() == mAllelePeptideData.size())
                        break;

                    continue;
                }

                BindData bindData = BindData.fromCsv(line, alleleIndex, peptideIndex, affinityIndex, predictedIndex, sourceIndex);

                if(bindData.Peptide.contains("X"))
                    continue;

                if(!allele.equals(currentAllele))
                {
                    currentAllele = allele;
                    currentBindList = Lists.newArrayList();
                    mAllelePeptideData.put(allele, currentBindList);
                }

                currentBindList.add(bindData);

                mDistinctPeptideLengths.add(bindData.peptideLength());
            }

            NE_LOGGER.info("loaded {} alleles with {} training data items from file({}) maxPeptideLength({})",
                    mAllelePeptideData.size(), mAllelePeptideData.values().stream().mapToInt(x -> x.size()).sum(),
                    mConfig.TrainingDataFile, getMaxPeptideLength());
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read training binding data file: {}", e.toString());
            return false;
        }

        return true;
    }

    private int getMaxPeptideLength() { return mDistinctPeptideLengths.stream().mapToInt(x -> x.intValue()).max().orElse(0); }

    private boolean loadRandomPredictionsData()
    {
        if(mConfig.RandomPeptidePredictionsFile == null)
            return true; // not required

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.RandomPeptidePredictionsFile));
            String header = fileReader.readLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

            int alleleIndex = fieldsIndexMap.get("Allele");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int predAffinityIndex = fieldsIndexMap.get("PredictedAffinity");

            String currentAllele = "";
            List<BindData> currentBindList = null;

            String line = "";

            int loadedAlleles = 0;
            int randomPeptides = 0;

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = items[alleleIndex];

                if(!mConfig.SpecificAlleles.isEmpty() && !mConfig.SpecificAlleles.contains(allele))
                {
                    if(mConfig.SpecificAlleles.size() == loadedAlleles)
                        break;

                    continue;
                }

                if(!mAllelePeptideData.containsKey(allele))
                    continue;

                BindData bindData = BindData.fromCsv(line, alleleIndex, peptideIndex, predAffinityIndex, mConfig.Constants.MaxAffinity);

                if(bindData.Peptide.contains("X"))
                    continue;

                if(!allele.equals(currentAllele))
                {
                    ++loadedAlleles;
                    currentAllele = allele;
                    currentBindList = mAllelePeptideData.get(allele);
                }

                currentBindList.add(bindData);
                ++randomPeptides;
            }

            NE_LOGGER.info("loaded {} random data items from file({})",
                    randomPeptides, mConfig.RandomPeptidePredictionsFile);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read random peptide binding data file: {}", e.toString());
            return false;
        }

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);
        options.addOption(AMINO_ACID_FREQ_FILE, true, "Amino acid frequency from proteome");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        BindTrainer neoBinder = new BindTrainer(cmd);
        neoBinder.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
