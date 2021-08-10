package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCountData.writeCounts;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.BindScoreMatrix.initMatrixWriter;
import static com.hartwig.hmftools.neo.bind.BindScoreMatrix.writeMatrixData;
import static com.hartwig.hmftools.neo.bind.GlobalWeights.GLOBAL_COUNTS;
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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class BindTrainer
{
    private final BinderConfig mConfig;

    private final Map<String,Map<Integer,List<BindData>>> mAllelePeptideData;

    private final Map<String,Map<Integer,BindCountData>> mAlleleBindCounts; // counts data by peptide length>
    private final Map<String,Map<Integer,BindScoreMatrix>> mAlleleBindMatrices;

    private final Set<Integer> mDistinctPeptideLengths;

    private final PosWeightModel mPosWeightModel;
    private final HlaSequences mHlaSequences;

    public BindTrainer(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mAllelePeptideData = Maps.newHashMap();
        mDistinctPeptideLengths = Sets.newHashSet();

        mHlaSequences = new HlaSequences();
        mHlaSequences.load(cmd.getOptionValue(HLA_DEFINITIONS_FILE));

        mPosWeightModel = new PosWeightModel(mConfig.Constants, mHlaSequences);

        mAlleleBindCounts = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
    }

    public void run()
    {
        NE_LOGGER.info("running NeoBinder on {} alleles", mConfig.RequiredAlleles.isEmpty() ? "all" : mConfig.RequiredAlleles.size());

        if(!loadTrainingData() || !loadRandomPredictionsData())
        {
            System.exit(1);
        }

        buildBindCountsData();
        buildPositionWeightMatrices();

        if(mConfig.RandomPeptides.WriteRandomDistribution || mConfig.RunScoring)
        {
            RandomPeptideDistribution randomDistribution = new RandomPeptideDistribution(mConfig.RandomPeptides);

            if(!randomDistribution.loadData())
            {
                randomDistribution.buildDistribution(mAlleleBindMatrices);
            }

            if(mConfig.RunScoring)
            {
                BindScorer scorer = new BindScorer(mConfig, mAllelePeptideData, mAlleleBindMatrices, randomDistribution);
                scorer.runScoring();

                if(mConfig.WriteLikelihood)
                {
                    BindingLikelihood compBinding = new BindingLikelihood();
                    final String relativeLikelihoodFile = mConfig.formFilename("rel_likelihood");
                    compBinding.buildAllelePeptideLikelihoods(mAllelePeptideData, relativeLikelihoodFile);
                }
            }
        }

        NE_LOGGER.info("NeoBinder complete");
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
        // fill in an gaps in allles or peptide lengths if they are required
        if(!mConfig.RequiredAlleles.isEmpty())
        {
            mConfig.RequiredAlleles.stream()
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

            for(BindCountData bindCounts : pepLenBindCounts)
            {
                totalAlelleCount += bindCounts.totalBindCount();
                mPosWeightModel.buildWeightedCounts(bindCounts, pepLenBindCounts);
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
                    mPosWeightModel.buildFinalWeightedCounts(bindCounts, allBindCounts, alleleTotalCounts);
                }
            }
        }

        BufferedWriter matrixWriter = mConfig.WritePosWeightMatrix || mConfig.WriteBindCounts ?
                initMatrixWriter(mConfig.formFilename("matrix_data"), getMaxPeptideLength(), mConfig.WriteBindCounts) : null;

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

        if(mPosWeightModel.getGlobalWeights().enabled())
        {
            mPosWeightModel.getGlobalWeights().writeGlobalCounts(
                    matrixWriter, getMaxPeptideLength(), mConfig.WritePosWeightMatrix, mConfig.WriteBindCounts);

            mAlleleBindMatrices.put(GLOBAL_COUNTS, mPosWeightModel.getGlobalWeights().getMatrixMap());
        }

        closeBufferedWriter(matrixWriter);
    }

    private boolean loadTrainingData()
    {
        if(!loadBindData(
                mConfig.TrainingDataFile, true, mConfig.RequiredAlleles, mConfig.RequiredPeptideLengths, mAllelePeptideData))
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
                mConfig.RandomPeptidePredictionsFile, false, mConfig.RequiredAlleles, mConfig.RequiredPeptideLengths, mAllelePeptideData))
        {
            return false;
        }

        if(mConfig.RandomPeptidePredictionsFile == null)
            return true; // not required

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);

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
