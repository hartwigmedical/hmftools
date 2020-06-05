package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class KnownFusionCache
{
    private final List<KnownFusionData> mData;
    private final Map<KnownFusionType,List<KnownFusionData>> mDataByType;

    // new-style combined input file
    public static final String KNOWN_FUSIONS_FILE = "known_fusion_file";

    // old-style file names
    public static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    public static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    public static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";

    private static final String FILE_DELIMITER = ",";

    private static final Logger LOGGER = LogManager.getLogger(KnownFusionCache.class);

    public KnownFusionCache()
    {
        mData = Lists.newArrayList();
        mDataByType = Maps.newHashMap();

        // initialise to avoid having to check for null
        Arrays.stream(KnownFusionType.values()).filter(x -> x != NONE).forEach(x -> mDataByType.put(x, Lists.newArrayList()));
    }

    public final List<KnownFusionData> getData() { return mData; }
    public final List<KnownFusionData> getDataByType(final KnownFusionType type) { return mDataByType.get(type); }

    public boolean hasKnownFusion(final String fiveGene, final String threeGene)
    {
        return mDataByType.get(KNOWN_PAIR).stream().anyMatch(x -> x.FiveGene.equals(fiveGene) && x.ThreeGene.equals(threeGene));
    }

    public boolean hasPromiscuousFiveGene(final String gene)
    {
        return mDataByType.get(PROMISCUOUS_5).stream().anyMatch(x -> x.FiveGene.equals(gene));
    }

    public boolean hasPromiscuousThreeGene(final String gene)
    {
        return mDataByType.get(PROMISCUOUS_3).stream().anyMatch(x -> x.ThreeGene.equals(gene));
    }

    public boolean isExonDelDup(final String geneName, final String transName, int fusedExonUp, int fusedExonDown)
    {
        for(final KnownFusionData knownData : mDataByType.get(EXON_DEL_DUP))
        {
            if(!knownData.FiveGene.equals(geneName) || !knownData.specificTransName().equals(transName))
                continue;

            if(fusedExonUp >= knownData.minFusedExons()[FS_UPSTREAM] && fusedExonUp <= knownData.maxFusedExons()[FS_UPSTREAM]
            && fusedExonDown >= knownData.minFusedExons()[FS_DOWNSTREAM] && fusedExonDown <= knownData.maxFusedExons()[FS_DOWNSTREAM])
            {
                return true;
            }
        }

        return false;
    }

    public boolean withinIgRegion(final String chromosome, int position)
    {
        if(mDataByType.get(IG_KNOWN_PAIR).stream().anyMatch(x -> x.withinIgRegion(chromosome, position)))
        {
            return true;
        }
        else if(mDataByType.get(IG_PROMISCUOUS).stream().anyMatch(x -> x.withinIgRegion(chromosome, position)))
        {
            return true;
        }

        return false;
    }

    public boolean matchesIgGene(final String chromosome, int position, byte orientation)
    {
        if(mDataByType.get(IG_KNOWN_PAIR).stream().anyMatch(x -> x.matchesIgGene(chromosome, position, orientation)))
        {
            return true;
        }
        else if(mDataByType.get(IG_PROMISCUOUS).stream().anyMatch(x -> x.matchesIgGene(chromosome, position, orientation)))
        {
            return true;
        }

        return false;
    }

    public boolean loadFromFile(@NotNull final CommandLine cmd)
    {
        if(cmd.hasOption(KNOWN_FUSIONS_FILE))
        {
            if(!loadFile(cmd.getOptionValue(KNOWN_FUSIONS_FILE)))
                return false;
        }
        else if(cmd.hasOption(PROMISCUOUS_THREE_CSV) || cmd.hasOption(PROMISCUOUS_FIVE_CSV) || cmd.hasOption(FUSION_PAIRS_CSV))
        {
            try
            {
                if(cmd.hasOption(PROMISCUOUS_THREE_CSV))
                {
                    final String threeGenesFile = cmd.getOptionValue(PROMISCUOUS_THREE_CSV);
                    loadOldStyleFiles(threeGenesFile, PROMISCUOUS_3);
                    LOGGER.info("loaded {} 3-prime promiscous genes from file: {}", mDataByType.get(PROMISCUOUS_3).size(), threeGenesFile);
                }

                if(cmd.hasOption(PROMISCUOUS_FIVE_CSV))
                {
                    final String fiveGenesFile = cmd.getOptionValue(PROMISCUOUS_FIVE_CSV);
                    loadOldStyleFiles(fiveGenesFile, PROMISCUOUS_5);
                    LOGGER.info("loaded {} 5-prime promiscous genes from file: {}", mDataByType.get(PROMISCUOUS_5).size(), fiveGenesFile);
                }

                if(cmd.hasOption(FUSION_PAIRS_CSV))
                {
                    final String knownPairsFile = cmd.getOptionValue(FUSION_PAIRS_CSV);
                    loadOldStyleFiles(knownPairsFile, KNOWN_PAIR);
                    LOGGER.info("loaded {} known fusion gene-pairs from file: {}", mDataByType.get(KNOWN_PAIR).size(), knownPairsFile);
                }
            }
            catch (IOException e)
            {
                LOGGER.warn("failed to load known fusion data: {}", e.toString());
                return false;
            }
        }

        for(Map.Entry<KnownFusionType,List<KnownFusionData>> entry : mDataByType.entrySet())
        {
            if(!entry.getValue().isEmpty())
            {
                LOGGER.info("loaded {} {} known-fusion records", entry.getValue().size(), entry.getKey());
            }
        }

        return true;
    }

    public void addData(final KnownFusionData data)
    {
        mData.add(data);
        mDataByType.get(data.Type).add(data);
    }

    private boolean loadFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LOGGER.error("file({}) not found", filename);
            return false;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return false;

            final Map<String,Integer> fieldIndexMap = createFieldsIndexMap(fileContents.get(0), FILE_DELIMITER);
            fileContents.remove(0);

            for(String data : fileContents)
            {
                KnownFusionData knownFusionData = KnownFusionData.fromCsv(data, fieldIndexMap);

                if(knownFusionData == null)
                {
                    LOGGER.error("file({}) invalid known fusion data will be skipped: {}", filename, data);
                    continue;
                }

                addData(knownFusionData);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("file({}) invalid known fusion data: {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private void loadOldStyleFiles(final String filename, final KnownFusionType type) throws IOException
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LOGGER.error("file({}) not found", filename);
            throw new IOException();
        }

        final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

        if(!fileContents.isEmpty())
        {
            // assumes a header row
            fileContents.remove(0);
        }

        for (String data : fileContents)
        {
            if(type == KNOWN_PAIR)
            {
                final String[] genes = data.split(FILE_DELIMITER);
                if(genes.length != 2)
                {
                    LOGGER.error("file({}) invalid known-pair data: {}", filename, data);
                    continue;
                }

                addData(new KnownFusionData(KNOWN_PAIR, genes[0], genes[1], "", "", ""));
            }
            else if(type == PROMISCUOUS_5)
            {
                addData(new KnownFusionData(PROMISCUOUS_5, data, "", "", "", ""));
            }
            else if(type == PROMISCUOUS_3)
            {
                addData(new KnownFusionData(PROMISCUOUS_3, "", data, "", "", ""));
            }
        }
    }

}
