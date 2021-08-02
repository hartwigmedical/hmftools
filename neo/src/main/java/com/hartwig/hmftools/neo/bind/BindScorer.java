package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

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

public class BindScorer
{
    private final BinderConfig mConfig;

    private final Map<String,List<BindData>> mAllelePeptideData;
    private final Map<String,Map<Integer,BindScoreMatrix>> mAlleleBindMatrices;

    public BindScorer(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);

        mAllelePeptideData = Maps.newHashMap();
        mAlleleBindMatrices = Maps.newHashMap();
    }

    public void run()
    {
        NE_LOGGER.info("running BindScorer");



        NE_LOGGER.info("BindScorer complete");
    }

    private void processingBindingData()
    {
    }

    private boolean loadBindData()
    {
        try
        {
            final List<String> lines = Files.readAllLines(new File(mConfig.TrainingDataFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get("Allele");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int typeIndex = fieldsIndexMap.get("AlleleType");

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

                BindData bindData = null; // BindData.fromCsv(line, alleleIndex, peptideIndex, affinityIndex, predictedIndex, sourceIndex);

                if(bindData.Peptide.contains("X"))
                    continue;

                if(!allele.equals(currentAllele))
                {
                    currentAllele = allele;
                    currentBindList = Lists.newArrayList();
                    mAllelePeptideData.put(allele, currentBindList);
                }

                currentBindList.add(bindData);
            }

            NE_LOGGER.info("loaded {} alleles with {} training data items from file({})",
                    mAllelePeptideData.size(), mAllelePeptideData.values().stream().mapToInt(x -> x.size()).sum(),
                    mConfig.TrainingDataFile);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read training binding data file: {}", e.toString());
            return false;
        }

        return true;
    }


    private void loadScoreMatrixData()
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        BindScorer bindScorer = new BindScorer(cmd);
        bindScorer.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
