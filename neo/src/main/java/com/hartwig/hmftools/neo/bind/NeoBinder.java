package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MAX_PEPTIDE_POSITIONS;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class NeoBinder
{
    private final BinderConfig mConfig;

    private final List<BindData> mBindDataList;

    public NeoBinder(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mBindDataList = Lists.newArrayList();
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
        BindMatrix simpleMatrix = new BindMatrix(MAX_PEPTIDE_POSITIONS);

        for(BindData bindData : mBindDataList)
        {
            double levelScore = deriveLevelScore(bindData.Affinity);
            simpleMatrix.processBindData(bindData, levelScore);
        }

        // write results
        String filename = mConfig.formFilename("freq_score");
        simpleMatrix.writeFrequencyData(filename);
    }

    private double deriveLevelScore(final double affinity)
    {
        for(double[] levelScore : mConfig.BindingLevelScores)
        {
            if(affinity < levelScore[0])
                return levelScore[1];
        }

        return mConfig.BindingLevelScores.get(mConfig.BindingLevelScores.size() - 1)[1];
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

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                BindData bindData = BindData.fromCsv(line, alleleIndex, peptideIndex, affinityIndex, otherInfoIndex);

                if(!mConfig.SpecificAlleles.isEmpty() && !mConfig.SpecificAlleles.contains(bindData.Allele))
                    continue;

                mBindDataList.add(bindData);
            }

            NE_LOGGER.info("loaded {} training data items from file({})", mBindDataList.size(), filename);
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
