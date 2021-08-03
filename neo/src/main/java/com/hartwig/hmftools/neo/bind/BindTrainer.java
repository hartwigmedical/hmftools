package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCountData.initMatrixWriter;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;
import static com.hartwig.hmftools.neo.bind.HlaSequences.POSITION_HLA_AA_FILE;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.LIKELY_INCORRECT;
import static com.hartwig.hmftools.neo.bind.PeptideWriteType.TRAINING;
import static com.hartwig.hmftools.neo.bind.BindCountData.initFrequencyWriter;
import static com.hartwig.hmftools.neo.utils.AminoAcidFrequency.AMINO_ACID_FREQ_FILE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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

    private final Map<String,Map<Integer,List<BindData>>> mAllelePeptideData;
    private final AminoAcidFrequency mAminoAcidFrequency;

    private final Map<String,Map<Integer,BindCountData>> mAlleleBindCounts; // counts data by peptide length>
    private final Map<String,Map<Integer,BindScoreMatrix>> mAlleleBindMatrices;

    private final Set<Integer> mDistinctPeptideLengths;

    private final BlosumMapping mBlosumMapping;
    private final HlaSequences mHlaSequences;

    public BindTrainer(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mAllelePeptideData = Maps.newHashMap();
        mDistinctPeptideLengths = Sets.newHashSet();

        mAminoAcidFrequency = new AminoAcidFrequency(cmd);
        mAminoAcidFrequency.loadFrequencies();

        mBlosumMapping = new BlosumMapping();

        mHlaSequences = new HlaSequences();
        mHlaSequences.load(cmd.getOptionValue(POSITION_HLA_AA_FILE), cmd.getOptionValue(HLA_DEFINITIONS_FILE));

        mAlleleBindCounts = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
    }

    public void run()
    {
        NE_LOGGER.info("running NeoBinder on {} alleles", mConfig.SpecificAlleles.isEmpty() ? "all" : mConfig.SpecificAlleles.size());

        if(!loadTrainingData() || !loadRandomPredictionsData())
        {
            System.exit(1);
        }

        if(mConfig.WritePosWeightMatrix && mAminoAcidFrequency.getAminoAcidFrequencies().isEmpty())
        {
            NE_LOGGER.warn("no amino acid frequencies loaded");
            System.exit(1);
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
            List<BindScoreMatrix> matrixList = BindScoreMatrix.loadFromCsv(mConfig.BindMatrixFile);

            for(BindScoreMatrix matrix : matrixList)
            {
                Map<Integer,BindScoreMatrix> alleleMap = mAlleleBindMatrices.get(matrix.Allele);

                if(alleleMap == null)
                {
                    alleleMap = Maps.newHashMap();
                    mAlleleBindMatrices.put(matrix.Allele, alleleMap);
                }

                alleleMap.put(matrix.PeptideLength, matrix);
            }
        }
        else
        {
            buildBindCountsData();
            buildPositionWeightMatrices();
        }
    }

    private void buildBindCountsData()
    {
        BufferedWriter freqWriter = mConfig.WriteFrequencyData ?
                initFrequencyWriter(mConfig.formFilename("pos_frequency")) : null;

        BufferedWriter pairWriter = mConfig.CalcPairs ?
                ComboCorrelations.initPairDataWriter(mConfig.formFilename("pair_score_prob")) : null;

        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();

            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer,BindCountData> pepLenCountsMap = Maps.newHashMap(); // counts data by peptide length
            mAlleleBindCounts.put(allele, pepLenCountsMap);

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<BindData> bindDataList = pepLenEntry.getValue();

                NE_LOGGER.debug("allele({}) length({}) processing {} data items", allele, pepLenEntry.getKey(), bindDataList.size());

                for(BindData bindData : bindDataList)
                {
                    int peptideLength = bindData.peptideLength();

                    BindCountData bindCounts = pepLenCountsMap.get(bindData.peptideLength());

                    if(bindCounts == null)
                    {
                        bindCounts = new BindCountData(allele, peptideLength);
                        pepLenCountsMap.put(peptideLength, bindCounts);
                    }

                    bindCounts.processBindData(bindData, mConfig.CalcPairs, mConfig.Constants);
                }
            }

            for(BindCountData bindCounts : pepLenCountsMap.values())
            {
                // write results
                bindCounts.logStats();

                if(mConfig.WriteFrequencyData)
                    bindCounts.writeFrequencyData(freqWriter);

                if(mConfig.CalcPairs)
                    ComboCorrelations.writePairData(pairWriter, bindCounts);
            }
        }

        closeBufferedWriter(freqWriter);
        closeBufferedWriter(pairWriter);
    }

    private void buildPositionWeightMatrices()
    {
        final Map<String,Integer> alleleTotalCounts = Maps.newHashMap();

        for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final List<BindCountData> pepLenBindCounts = alleleEntry.getValue().values().stream().collect(Collectors.toList());
            int totalAlelleCount = 0;

            NE_LOGGER.debug("allele({}) build matrix for {} peptide lengths", allele, pepLenBindCounts.size());

            for(BindCountData bindCounts : pepLenBindCounts)
            {
                totalAlelleCount += bindCounts.totalBindCount();
                bindCounts.buildWeightedCounts(pepLenBindCounts, mConfig.Constants);
            }

            alleleTotalCounts.put(allele, totalAlelleCount);
        }

        Map<Integer,List<BindCountData>> countsByLength = Maps.newHashMap();
        mDistinctPeptideLengths.forEach(x -> countsByLength.put(x, Lists.newArrayList()));

        // now factor each allele's weighted counts into all the others
        for(Integer peptideLength : mDistinctPeptideLengths)
        {
            List<BindCountData> allBindCounts = Lists.newArrayList();
            countsByLength.put(peptideLength, allBindCounts);

            for(Map<Integer,BindCountData> alleleEntry : mAlleleBindCounts.values())
            {
                BindCountData bindCounts = alleleEntry.get(peptideLength);

                if(bindCounts != null)
                    allBindCounts.add(bindCounts);
            }

            // now apply the list of counts across all alleles per peptide length
            for(Map<Integer,BindCountData> alleleEntry : mAlleleBindCounts.values())
            {
                BindCountData bindCounts = alleleEntry.get(peptideLength);

                if(bindCounts != null)
                {
                    bindCounts.buildFinalWeightedCounts(allBindCounts, alleleTotalCounts, mConfig.Constants, mBlosumMapping, mHlaSequences);
                }
            }
        }

        BufferedWriter matrixWriter = mConfig.WritePosWeightMatrix ?
                initMatrixWriter(mConfig.formFilename("matrix_data"), getMaxPeptideLength(), mConfig.WriteBindCounts) : null;

        for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,BindCountData> peptideLengthCounts = alleleEntry.getValue();

            // calculate peptide length frequency per allele
            final Map<Integer,Integer> peptideLengthFrequency = Maps.newHashMap();
            peptideLengthCounts.entrySet().forEach(x -> peptideLengthFrequency.put(x.getKey(), x.getValue().totalBindCount()));

            Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = Maps.newHashMap();
            mAlleleBindMatrices.put(allele, peptideLengthMatrixMap);

            for(BindCountData bindCounts : peptideLengthCounts.values())
            {
                BindScoreMatrix matrix = bindCounts.createMatrix(mAminoAcidFrequency, peptideLengthFrequency);
                peptideLengthMatrixMap.put(matrix.PeptideLength, matrix);

                if(mConfig.WritePosWeightMatrix)
                    bindCounts.writeMatrixData(matrixWriter, matrix, getMaxPeptideLength(), mConfig.WriteBindCounts);
            }
        }

        closeBufferedWriter(matrixWriter);
    }

    private void runScoring()
    {
        RandomPeptideDistribution randomDistribution = new RandomPeptideDistribution(mConfig);

        if(!randomDistribution.hasData())
        {
            randomDistribution.buildDistribution(mAlleleBindMatrices);
        }

        NE_LOGGER.info("writing allele summaries");

        BufferedWriter alleleWriter = initAlleleSummaryWriter();

        // rank both the training and random peptide data using the newly created data and the random peptide distributions
        for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = mAlleleBindMatrices.get(allele);

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                final List<BindData> bindDataList = pepLenEntry.getValue();
                for(BindData bindData : bindDataList)
                {
                    BindScoreMatrix matrix = peptideLengthMatrixMap.get(bindData.peptideLength());

                    if(matrix == null)
                        continue;

                    double score = matrix.calcScore(bindData.Peptide);
                    bindData.setScoreData(score, randomDistribution.getScoreRank(allele, bindData.peptideLength(), score));
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
            writer.write("Allele,PeptideLength,Source,TrainingCount,RandomCount,AUC");
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
            final Map<Integer,List<BindData>> pepLenBindDataMap = mAllelePeptideData.get(allele);

            int trainingCount = 0;
            int randomCount = 0;
            int errors = 0;

            // internal model assesses AUC per peptide length
            // external models across an entire allele

            List<AucData> alleleAucMcfData = Lists.newArrayList();

            for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<AucData> alleleAucData = Lists.newArrayList();
                List<AucData> alleleAucRankData = Lists.newArrayList();

                for(BindData bindData : pepLenEntry.getValue())
                {
                    if(bindData.isTraining())
                        ++trainingCount;
                    else
                        ++randomCount;

                    boolean isPositive = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;

                    alleleAucData.add(new AucData(isPositive, bindData.score(), false));
                    alleleAucRankData.add(new AucData(isPositive, bindData.rankPercentile(), true));
                    alleleAucMcfData.add(new AucData(isPositive, bindData.affinityPercentile(), true));
                }

                // double auc = AucCalc.calcScoresAuc(alleleAucData, Level.TRACE);
                double aucPerc = AucCalc.calcPercentilesAuc(alleleAucRankData, Level.TRACE);

                NE_LOGGER.info(String.format("allele(%s) peptideLength(%d) peptides(train=%d, rand=%d) AUC(%.4f)",
                        allele, pepLenEntry.getKey(), trainingCount, randomCount, aucPerc));

                writer.write(String.format("%s,%d,%s,%d,%d,%.4f",
                        allele, pepLenEntry.getKey(), "MODEL", trainingCount, randomCount, aucPerc));
                writer.newLine();
            }

            double aucMcf = AucCalc.calcPercentilesAuc(alleleAucMcfData, Level.TRACE);

            NE_LOGGER.info(String.format("allele(%s) McFlurry peptides(train=%d, rand=%d) AUC(%.4f)",
                    allele, trainingCount, randomCount, aucMcf));

            writer.write(String.format("%s,%d,%s,%d,%d,%.4f",
                    allele, 0, "MCF", trainingCount, randomCount, aucMcf));
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

            for(Map.Entry<String,Map<Integer,List<BindData>>> alleleEntry : mAllelePeptideData.entrySet())
            {
                final String allele = alleleEntry.getKey();
                final Map<Integer,List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

                for(List<BindData> bindDataList : pepLenBindDataMap.values())
                {
                    for(BindData bindData : bindDataList)
                    {
                        if(mConfig.WritePeptideType == LIKELY_INCORRECT)
                        {
                            boolean expectBinds = bindData.Affinity < mConfig.Constants.BindingAffinityHigh;
                            boolean inTopPerc = bindData.rankPercentile() < 0.02;
                            if(expectBinds == inTopPerc)
                                continue;
                        }
                        else if(mConfig.WritePeptideType == TRAINING && !bindData.isTraining())
                        {
                            continue;
                        }

                        writer.write(String.format("%s,%s,%s,%.4f,%.4f,%.1f,%.1f",
                                allele, bindData.Peptide, bindData.Source,
                                bindData.score(), bindData.rankPercentile(), bindData.Affinity, bindData.predictedAffinity()));
                        writer.newLine();
                    }
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
        if(!loadBindData(
                mConfig.TrainingDataFile, true, mConfig.SpecificAlleles, mConfig.SpecificPeptideLengths, mAllelePeptideData))
        {
            return false;
        }

        for(Map<Integer,List<BindData>> pepLenBindDataMap : mAllelePeptideData.values())
        {
            pepLenBindDataMap.keySet().forEach(x -> mDistinctPeptideLengths.add(x));
        }

        NE_LOGGER.info("{} peptide lengths, maxPeptideLength({})",
                mDistinctPeptideLengths.size(), getMaxPeptideLength());

        return true;
    }

    private int getMaxPeptideLength() { return mDistinctPeptideLengths.stream().mapToInt(x -> x.intValue()).max().orElse(0); }

    private boolean loadRandomPredictionsData()
    {
        if(!loadBindData(
                mConfig.RandomPeptidePredictionsFile, false, mConfig.SpecificAlleles, mConfig.SpecificPeptideLengths, mAllelePeptideData))
        {
            return false;
        }

        if(mConfig.RandomPeptidePredictionsFile == null)
            return true; // not required

        /*
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
         */

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
