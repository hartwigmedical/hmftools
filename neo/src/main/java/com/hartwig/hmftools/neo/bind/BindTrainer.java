package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCountData.writeCounts;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.BindScoreMatrix.initMatrixWriter;
import static com.hartwig.hmftools.neo.bind.BindScoreMatrix.writeMatrixData;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_FLANK_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;
import static com.hartwig.hmftools.neo.bind.BindCountData.initFrequencyWriter;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.MatrixUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class BindTrainer
{
    private final TrainConfig mConfig;
    private final ScoreConfig mScoreConfig;

    private final Map<String,Map<Integer,List<BindData>>> mAllelePeptideData;

    private final Map<String,Map<Integer,BindCountData>> mAlleleBindCounts; // counts data by peptide length>
    private final Map<String,Map<Integer,BindScoreMatrix>> mAlleleBindMatrices;
    private final FlankCounts mFlankCounts;
    private final FlankScores mFlankScores;
    private final RecognitionCounts mRecognitionCounts;

    private final Set<Integer> mDistinctPeptideLengths;

    private final PosWeightModel mPosWeightModel;
    private final HlaSequences mHlaSequences;

    public BindTrainer(final CommandLine cmd)
    {
        mConfig = new TrainConfig(cmd);
        mScoreConfig = new ScoreConfig(cmd);
        mAllelePeptideData = Maps.newHashMap();
        mDistinctPeptideLengths = Sets.newHashSet();

        mHlaSequences = new HlaSequences();
        mHlaSequences.load(cmd.getOptionValue(HLA_DEFINITIONS_FILE));

        mPosWeightModel = new PosWeightModel(mConfig.Constants, mHlaSequences);

        mAlleleBindCounts = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
        mFlankCounts = new FlankCounts();
        mFlankScores = new FlankScores();
        mRecognitionCounts = new RecognitionCounts();
    }

    public void run()
    {
        NE_LOGGER.info("building training data from {} alleles for {} required alleles",
                mAllelePeptideData.size(), mConfig.RequiredOutputAlleles.size());

        if(!loadTrainingData())
        {
            System.exit(1);
        }

        // check that the HLA defintions cover the training and required (if specified) alleles if specified
        Set<String> allAlleles = mAllelePeptideData.keySet().stream().collect(Collectors.toSet());
        mConfig.RequiredOutputAlleles.forEach(x -> allAlleles.add(x));

        List<String> missingAlleleDefinitions = allAlleles.stream()
                .filter(x -> !mHlaSequences.hasAlleleDefinition(x)).collect(Collectors.toList());

        if(!missingAlleleDefinitions.isEmpty())
        {
            NE_LOGGER.error("HLA definitions are missing the following alleles: {}", missingAlleleDefinitions);
            System.exit(1);
        }

        buildBindCountsData();
        buildPositionWeightMatrices();

        if(mConfig.RandomPeptides.WriteRandomDistribution)
        {
            RandomPeptideDistribution randomDistribution = new RandomPeptideDistribution(mConfig.RandomPeptides);

            if(!randomDistribution.loadData())
                randomDistribution.buildDistribution(mAlleleBindMatrices, mFlankScores);

            if(mConfig.WriteLikelihood)
            {
                BindScorer scorer = new BindScorer(mScoreConfig, mAllelePeptideData, mAlleleBindMatrices, randomDistribution, mFlankScores);
                scorer.runScoring();

                BindingLikelihood bindingLikelihood = new BindingLikelihood();
                final String relativeLikelihoodFile = mConfig.formOutputFilename(FILE_ID_LIKELIHOOD);
                bindingLikelihood.buildAllelePeptideLikelihoods(mAllelePeptideData, relativeLikelihoodFile);

                if(mConfig.WritePanLengthDistribution)
                {
                    randomDistribution.buildLikelihoodDistribution(mAlleleBindMatrices, mFlankScores, bindingLikelihood);
                }
            }
        }

        NE_LOGGER.info("NeoBinder complete");
    }

    private void buildBindCountsData()
    {
        BufferedWriter freqWriter = mConfig.WriteFrequencyData ?
                initFrequencyWriter(mConfig.formOutputFilename("pos_frequency")) : null;

        BufferedWriter pairWriter = mConfig.CalcPairs ?
                ComboCorrelations.initPairDataWriter(mConfig.formOutputFilename("pair_score_prob")) : null;

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

                    bindCounts.processBindData(bindData, mConfig.CalcPairs);

                    if(mConfig.ApplyFlanks)
                        mFlankCounts.processBindData(bindData);
                }
            }

            for(BindCountData bindCounts : pepLenCountsMap.values())
            {
                // write results
                // bindCounts.logStats();

                if(mConfig.WriteFrequencyData)
                    bindCounts.writeFrequencyData(freqWriter);

                if(mConfig.CalcPairs)
                    ComboCorrelations.writePairData(pairWriter, bindCounts);
            }
        }

        if(mConfig.ApplyFlanks)
        {
            mFlankScores.createMatrix(mFlankCounts.getBindCounts());
        }

        if(mConfig.ApplyFlanks && (mConfig.WriteBindCounts || mConfig.WritePosWeightMatrix))
        {
            mFlankCounts.logStats();
            mFlankCounts.writeData(mConfig.formOutputFilename(FILE_ID_FLANK_POS_WEIGHT), mFlankScores);
        }

        closeBufferedWriter(freqWriter);
        closeBufferedWriter(pairWriter);
    }

    private void buildPositionWeightMatrices()
    {
        // fill in an gaps in alleles or peptide lengths if they are required for the training data output
        if(!mConfig.RequiredOutputAlleles.isEmpty())
        {
            mConfig.RequiredOutputAlleles.stream()
                    .filter(x -> !mAlleleBindCounts.containsKey(x))
                    .forEach(x -> mAlleleBindCounts.put(x, Maps.newHashMap()));
        }

        final Map<String,Integer> alleleTotalCounts = Maps.newHashMap();

        for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer,BindCountData> pepLenBindCountsMap = alleleEntry.getValue();

            if(!mConfig.RequiredPeptideLengths.isEmpty())
            {
                mConfig.RequiredPeptideLengths.stream()
                        .filter(x -> !pepLenBindCountsMap.containsKey(x))
                        .forEach(x -> pepLenBindCountsMap.put(x, new BindCountData(allele, x)));
            }

            final List<BindCountData> pepLenBindCounts = pepLenBindCountsMap.values().stream().collect(Collectors.toList());

            int totalAlelleCount = 0;

            NE_LOGGER.debug("allele({}) building counts data for {} peptide lengths", allele, pepLenBindCounts.size());

            if(mPosWeightModel.noiseEnabled())
            {
                for(BindCountData bindCounts : pepLenBindCounts)
                {
                    mPosWeightModel.buildNoiseCounts(bindCounts);
                }
            }

            // calculate weighted counts across peptide lengths per position and amino acid
            for(BindCountData bindCounts : pepLenBindCounts)
            {
                totalAlelleCount += bindCounts.totalBindCount();
                mPosWeightModel.buildWeightedCounts(bindCounts, pepLenBindCounts);
            }

            // calculate position totals across the lengths for use when blending alleles
            mPosWeightModel.buildPositionAdjustedTotals(pepLenBindCounts);

            alleleTotalCounts.put(allele, totalAlelleCount);
        }

        // now factor each allele's weighted counts into all the others
        for(Integer peptideLength : mDistinctPeptideLengths)
        {
            List<BindCountData> allBindCounts = Lists.newArrayList();

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
                    mPosWeightModel.buildFinalWeightedCounts(bindCounts, allBindCounts);
                }
            }
        }

        logCalcAlleles();

        BufferedWriter matrixWriter = mConfig.WritePosWeightMatrix || mConfig.WriteBindCounts ?
                initMatrixWriter(mConfig.formOutputFilename(FILE_ID_POS_WEIGHT), getMaxPeptideLength()) : null;

        for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();

            NE_LOGGER.debug("allele({}) creating matrix data", allele);

            final Map<Integer,BindCountData> pepLenBindCounts = alleleEntry.getValue();

            Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = Maps.newHashMap();
            mAlleleBindMatrices.put(allele, peptideLengthMatrixMap);

            for(BindCountData bindCounts : pepLenBindCounts.values())
            {
                BindScoreMatrix matrix = mPosWeightModel.createMatrix(bindCounts);
                peptideLengthMatrixMap.put(matrix.PeptideLength, matrix);

                if(mConfig.WritePosWeightMatrix)
                    writeMatrixData(matrixWriter, matrix, getMaxPeptideLength());

                if(mConfig.WriteBindCounts)
                    writeCounts(matrixWriter, bindCounts, getMaxPeptideLength(), mPosWeightModel.noiseEnabled());
            }
        }

        closeBufferedWriter(matrixWriter);
    }

    private void logCalcAlleles()
    {
        if(mConfig.LogCalcAlleles.isEmpty())
            return;

        BufferedWriter writer = initMatrixWriter(mConfig.formOutputFilename("allele_motif_calc"), getMaxPeptideLength());

        for(String allele : mConfig.LogCalcAlleles)
        {
            Map<Integer,BindCountData> pepLenBindCounts = mAlleleBindCounts.get(allele);

            for(Map.Entry<Integer,BindCountData> pepLenEntry : pepLenBindCounts.entrySet())
            {
                final BindCountData bindCounts = pepLenEntry.getValue();

                for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
                {
                    final String otherAllele = alleleEntry.getKey();

                    if(otherAllele.equals(allele))
                        continue;

                    final BindCountData otherBindCounts = alleleEntry.getValue().get(bindCounts.PeptideLength);

                    BindCountData bindCountsTemp = new BindCountData(bindCounts.Allele, bindCounts.PeptideLength);

                    MatrixUtils.copy(bindCounts.getBindCounts(), bindCountsTemp.getBindCounts());
                    MatrixUtils.copy(bindCounts.getNoiseCounts(), bindCountsTemp.getNoiseCounts());
                    MatrixUtils.copy(bindCounts.getWeightedCounts(), bindCountsTemp.getWeightedCounts());

                    mPosWeightModel.buildFinalWeightedCounts(bindCountsTemp, Lists.newArrayList(otherBindCounts));

                    BindCountData bindCountsLog = new BindCountData(
                            String.format("%s->%s", otherAllele, bindCounts.Allele), bindCounts.PeptideLength);

                    writeCounts(writer, bindCountsLog, "AlleleMotif", bindCountsTemp.getFinalWeightedCounts(), getMaxPeptideLength());
                }
            }
        }

        closeBufferedWriter(writer);
    }

    private boolean loadTrainingData()
    {
        if(!loadBindData(mConfig.TrainingDataFile, mConfig.RequiredPeptideLengths, mAllelePeptideData))
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        TrainConfig.addCmdLineArgs(options);

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
