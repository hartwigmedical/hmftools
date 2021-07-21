package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MAX_PEPTIDE_POSITIONS;
import static com.hartwig.hmftools.neo.bind.BindMatrix.initCombProbabilityWriter;
import static com.hartwig.hmftools.neo.bind.BindMatrix.initFrequencyWriter;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class NeoBinder
{
    private final BinderConfig mConfig;

    private final Map<String,List<BindData>> mAlleleBindData;

    private final int mBindAffinityThreshold;

    public NeoBinder(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mAlleleBindData = Maps.newHashMap();

        mBindAffinityThreshold = 500;
    }

    public void run()
    {
        if(!loadTrainingData(mConfig.TrainingDataFile))
        {
            System.exit(1);
        }

        processingBindingData();

        NE_LOGGER.info("NeoBinder complete");
    }

    private void processingBindingData()
    {
        BufferedWriter freqWriter = initFrequencyWriter(mConfig.formFilename("freq_score"));
        BufferedWriter probWriter = initCombProbabilityWriter(mConfig.formFilename("pois_prob"));

        for(Map.Entry<String,List<BindData>> entry : mAlleleBindData.entrySet())
        {
            String allele = entry.getKey();
            BindMatrix simpleMatrix = new BindMatrix(allele, MAX_PEPTIDE_POSITIONS);

            for(BindData bindData : entry.getValue())
            {
                double levelScore = deriveLevelScore(bindData.Affinity);
                simpleMatrix.processBindData(bindData, levelScore, bindData.Affinity < mBindAffinityThreshold);
            }

            // write results

            simpleMatrix.logStats();
            simpleMatrix.writeFrequencyData(allele, freqWriter, probWriter);
        }

        closeBufferedWriter(freqWriter);
        closeBufferedWriter(probWriter);
    }

    private double deriveLevelScore(final double affinity)
    {
        if(!mConfig.BindingLevelScores.isEmpty())
        {
            for(double[] levelScore : mConfig.BindingLevelScores)
            {
                if(affinity < levelScore[0])
                    return levelScore[1];
            }

            return mConfig.BindingLevelScores.get(mConfig.BindingLevelScores.size() - 1)[1];
        }

        if(affinity >= mConfig.BindingLevelExponent)
            return 0;

        if(affinity <= 0)
            return 1;

        return 1 - log(mConfig.BindingLevelExponent, affinity);
    }

    private boolean loadTrainingData(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get("Allele");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int affinityIndex = fieldsIndexMap.get("Affinity");
            int otherInfoIndex = fieldsIndexMap.get("OtherInfo");

            String currentAllele = "";
            List<BindData> currentBindList = null;

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                String allele = items[alleleIndex];

                if(!mConfig.SpecificAlleles.isEmpty() && !mConfig.SpecificAlleles.contains(allele))
                {
                    if(mConfig.SpecificAlleles.size() == mAlleleBindData.size())
                        break;

                    continue;
                }

                BindData bindData = BindData.fromCsv(line, alleleIndex, peptideIndex, affinityIndex, otherInfoIndex);

                if(!allele.equals(currentAllele))
                {
                    currentAllele = allele;
                    currentBindList = Lists.newArrayList();
                    mAlleleBindData.put(allele, currentBindList);
                }

                currentBindList.add(bindData);
            }

            NE_LOGGER.info("loaded {} alleles with {} training data items from file({})",
                    mAlleleBindData.size(), mAlleleBindData.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read training binding data file: {}", e.toString());
            return false;
        }

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        NeoBinder neoBinder = new NeoBinder(cmd);
        neoBinder.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
