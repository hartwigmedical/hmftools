package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR_UNMAPPABLE_3;
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
import java.util.StringJoiner;

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

    private final List<KnownFusionData> mIgRegionData;

    public static final String KNOWN_FUSIONS_FILE = "known_fusion_file";
    private static final String FILE_DELIMITER = ",";

    private static final Logger LOGGER = LogManager.getLogger(KnownFusionCache.class);

    public KnownFusionCache()
    {
        mData = Lists.newArrayList();
        mDataByType = Maps.newHashMap();
        mIgRegionData = Lists.newArrayList();

        // initialise to avoid having to check for null
        Arrays.stream(KnownFusionType.values()).filter(x -> x != NONE).forEach(x -> mDataByType.put(x, Lists.newArrayList()));
    }

    public final List<KnownFusionData> getData() { return mData; }
    public final List<KnownFusionData> getDataByType(final KnownFusionType type) { return mDataByType.get(type); }

    public boolean hasKnownFusion(final String fiveGene, final String threeGene)
    {
        return mDataByType.get(KNOWN_PAIR).stream().anyMatch(x -> x.FiveGene.equals(fiveGene) && x.ThreeGene.equals(threeGene));
    }

    public boolean hasKnownUnmappable3Fusion(final String fiveGene, final String threeGene)
    {
        return mDataByType.get(KNOWN_PAIR_UNMAPPABLE_3).stream().anyMatch(x -> x.FiveGene.equals(fiveGene) && x.ThreeGene.equals(threeGene));
    }

    public boolean matchesKnownFusionGene(final GeneAnnotation gene)
    {
        if(mDataByType.get(KNOWN_PAIR).stream()
            .anyMatch(x -> (gene.isUpstream() && x.FiveGene.equals(gene.GeneName))
                    || (!gene.isUpstream() && x.ThreeGene.equals(gene.GeneName))))
        {
            return true;
        }

        if(mDataByType.get(KNOWN_PAIR_UNMAPPABLE_3).stream()
                .anyMatch(x -> !gene.isUpstream() && x.ThreeGene.equals(gene.GeneName)))
        {
            return true;
        }

        return false;
    }

    public boolean hasPromiscuousFiveGene(final String gene)
    {
        return mDataByType.get(PROMISCUOUS_5).stream().anyMatch(x -> x.FiveGene.equals(gene));
    }

    public boolean hasPromiscuousThreeGene(final String gene)
    {
        return mDataByType.get(PROMISCUOUS_3).stream().anyMatch(x -> x.ThreeGene.equals(gene));
    }

    public boolean isExonDelDupTrans(final String transName)
    {
        return mDataByType.get(EXON_DEL_DUP).stream().anyMatch(x -> x.specificTransName().equals(transName));
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
        return mIgRegionData.stream().anyMatch(x -> x.withinIgRegion(chromosome, position));
    }

    public boolean matchesIgGene(final String chromosome, int position, byte orientation)
    {
        return mIgRegionData.stream().anyMatch(x -> x.matchesIgGene(chromosome, position, orientation));
    }

    public boolean loadFromFile(@NotNull final CommandLine cmd)
    {
        if(!cmd.hasOption(KNOWN_FUSIONS_FILE))
            return true;

        if(!loadFile(cmd.getOptionValue(KNOWN_FUSIONS_FILE)))
            return false;

        StringJoiner refDataStr = new StringJoiner(", ");

        for(Map.Entry<KnownFusionType,List<KnownFusionData>> entry : mDataByType.entrySet())
        {
            if(!entry.getValue().isEmpty())
            {
                refDataStr.add(String.format("%s(%d)", entry.getKey(), entry.getValue().size()));
            }
        }

        LOGGER.info("loaded known fusion data: {}", refDataStr.toString());
        return true;
    }

    public void addData(final KnownFusionData data)
    {
        mData.add(data);
        mDataByType.get(data.Type).add(data);

        if(data.igRegion() != null)
            mIgRegionData.add(data);
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
                try
                {
                    KnownFusionData knownFusionData = KnownFusionData.fromCsv(data, fieldIndexMap);
                    addData(knownFusionData);
                }
                catch (Exception e)
                {
                    LOGGER.error("file({}) invalid known fusion data will be skipped: {}", filename, data);
                }
            }
        }
        catch (IOException e)
        {
            LOGGER.error("file({}) invalid known fusion data: {}", filename, e.toString());
            return false;
        }

        return true;
    }
}
