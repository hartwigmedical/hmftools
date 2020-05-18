package com.hartwig.hmftools.common.fusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class KnownFusionData
{

    private final List<String[]> mKnownPairs;
    private final List<String> mPromiscuousThreeGenes;
    private final List<String> mPromiscuousFiveGenes;

    public static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    public static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    public static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";

    public static final int FIVE_GENE = 0;
    public static final int THREE_GENE = 1;
    private static final String FILE_DELIMITER = ",";

    private static final Logger LOGGER = LogManager.getLogger(KnownFusionData.class);

    public KnownFusionData()
    {
        mKnownPairs = Lists.newArrayList();
        mPromiscuousThreeGenes = Lists.newArrayList();
        mPromiscuousFiveGenes = Lists.newArrayList();
    }

    public final List<String[]> knownPairs() { return mKnownPairs; }
    public final List<String> promiscuousThreeGenes() { return mPromiscuousThreeGenes; }
    public final List<String> promiscuousFiveGenes() { return mPromiscuousFiveGenes; }

    public boolean hasKnownFusion(final String fiveGene, final String threeGene)
    {
        return mKnownPairs.stream().anyMatch(x -> x[FIVE_GENE].equals(fiveGene) && x[THREE_GENE].equals(threeGene));
    }

    public boolean hasPromiscuousFiveGene(final String gene) { return mPromiscuousFiveGenes.contains(gene); }
    public boolean hasPromiscuousThreeGene(final String gene) { return mPromiscuousThreeGenes.contains(gene); }

    public boolean intergenicPromiscuousMatch(final String fiveGene, final String threeGene)
    {
        // cannot be same gene
        if(fiveGene.equals(threeGene))
            return false;

        return hasPromiscuousFiveGene(fiveGene) || hasPromiscuousThreeGene(threeGene);
    }

    public boolean intragenicPromiscuousMatch(final String fiveGene, final String threeGene)
    {
        return fiveGene.equals(threeGene) && hasPromiscuousThreeGene(threeGene);
    }

    public boolean loadFromFile(@NotNull final CommandLine cmd)
    {
        if(!cmd.hasOption(PROMISCUOUS_THREE_CSV) && !cmd.hasOption(PROMISCUOUS_FIVE_CSV) && !cmd.hasOption(FUSION_PAIRS_CSV))
            return false;

        try
        {
            if(cmd.hasOption(PROMISCUOUS_THREE_CSV))
            {
                final String threeGenesFile = cmd.getOptionValue(PROMISCUOUS_THREE_CSV);
                mPromiscuousThreeGenes.addAll(loadFile(threeGenesFile, 1));
                LOGGER.info("loaded {} 3-prime promiscous genes from file: {}", mPromiscuousThreeGenes.size(), threeGenesFile);
            }

            if(cmd.hasOption(PROMISCUOUS_FIVE_CSV))
            {
                final String fiveGenesFile = cmd.getOptionValue(PROMISCUOUS_FIVE_CSV);
                mPromiscuousFiveGenes.addAll(loadFile(fiveGenesFile, 1));
                LOGGER.info("loaded {} 5-prime promiscous genes from file: {}", mPromiscuousFiveGenes.size(), fiveGenesFile);
            }

            if(cmd.hasOption(FUSION_PAIRS_CSV))
            {
                final String knownPairsFile = cmd.getOptionValue(FUSION_PAIRS_CSV);
                final List<String> knownPairs = loadFile(knownPairsFile, 2);
                mKnownPairs.addAll(knownPairs.stream().map(x -> x.split(FILE_DELIMITER)).collect(Collectors.toList()));
                LOGGER.info("loaded {} known fusion gene-pairs from file: {}", mKnownPairs.size(), knownPairsFile);
            }
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load known fusion data: {}", e.toString());
            return false;
        }

        return true;
    }

    private final List<String> loadFile(final String filename, int expectedDataCount) throws IOException
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LOGGER.error("file({}) not found", filename);
            throw new IOException();
        }

        List<String> fileContents = Files.readAllLines(new File(filename).toPath());

        if(!fileContents.isEmpty())
        {
            // assumes a header row
            fileContents.remove(0);
        }

        if(!fileContents.isEmpty())
        {
            int index = 0;

            while (index < fileContents.size())
            {
                if (fileContents.get(index).split(FILE_DELIMITER).length != expectedDataCount)
                {
                    LOGGER.error("file({}) invalid data will be skipped: {}", filename, fileContents.get(index));
                    fileContents.remove(index);
                }
                else
                {
                    ++index;
                }
            }
        }

        return fileContents;
    }

}
