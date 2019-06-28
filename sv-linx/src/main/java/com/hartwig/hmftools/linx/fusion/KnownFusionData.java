package com.hartwig.hmftools.linx.fusion;

import static java.nio.file.Files.readAllLines;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;

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
        for(String[] pair : mKnownPairs)
        {
            if(pair[FIVE_GENE].equals(fiveGene) && pair[THREE_GENE].equals(threeGene))
                return true;
        }

        return false;
    }

    public boolean hasPromiscuousFiveGene(final String gene) { return mPromiscuousFiveGenes.contains(gene); }
    public boolean hasPromiscuousThreeGene(final String gene) { return mPromiscuousThreeGenes.contains(gene); }

    public boolean intergenicPromiscuousMatch(final String fiveGene, final String threeGene)
    {
        // cannot be same gene
        if(fiveGene.equals(threeGene))
            return false;

        return hasPromiscuousFiveGene(fiveGene) || hasPromiscuousThreeGene(threeGene);

        // return fiveGeneNames.stream().noneMatch(threeGeneNames::contains) && (
        //  fiveGeneNames.stream().anyMatch(promiscuousFive()::containsKey) || threeGeneNames.stream()
        //  .anyMatch(promiscuousThree()::containsKey));
    }

    public boolean intragenicPromiscuousMatch(final String fiveGene, final String threeGene)
    {
        return fiveGene.equals(threeGene) && hasPromiscuousThreeGene(threeGene);

        // return fiveGeneNames.stream().anyMatch(threeGeneNames::contains) && threeGeneNames.stream()
        //  .anyMatch(promiscuousThree()::containsKey);
    }

    public boolean loadFromFile(@NotNull final CommandLine cmd)
    {
        try
        {
            if(cmd.hasOption(PROMISCUOUS_THREE_CSV))
            {
                final String threeGenesFile = cmd.getOptionValue(PROMISCUOUS_THREE_CSV);
                mPromiscuousThreeGenes.addAll(loadFile(threeGenesFile, 1));
            }

            if(cmd.hasOption(PROMISCUOUS_FIVE_CSV))
            {
                final String fiveGenesFile = cmd.getOptionValue(PROMISCUOUS_FIVE_CSV);
                mPromiscuousFiveGenes.addAll(loadFile(fiveGenesFile, 1));
            }

            if(cmd.hasOption(FUSION_PAIRS_CSV))
            {
                final String knownPairsFile = cmd.getOptionValue(FUSION_PAIRS_CSV);
                final List<String> knownPairs = loadFile(knownPairsFile, 2);
                mKnownPairs.addAll(knownPairs.stream().map(x -> x.split(FILE_DELIMITER)).collect(Collectors.toList()));
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
            return Lists.newArrayList();

        List<String> fileContents = Files.readAllLines(new File(filename).toPath());

        if(!fileContents.isEmpty())
        {
            // remove header
            fileContents.remove(0);
        }

        if(!fileContents.isEmpty())
        {
            if(fileContents.get(0).split(FILE_DELIMITER).length > expectedDataCount)
            {
                for (int i = 0; i < fileContents.size(); ++i)
                {
                    String[] items = fileContents.get(i).split(FILE_DELIMITER);
                    String trimmedStr = items[0].replaceAll("\"", "");

                    for (int j = 1; j < expectedDataCount; ++j)
                    {
                        trimmedStr += FILE_DELIMITER + items[j].replaceAll("\"", "");
                    }

                    fileContents.set(i, trimmedStr);
                }
            }
        }

        return fileContents;
    }

}
