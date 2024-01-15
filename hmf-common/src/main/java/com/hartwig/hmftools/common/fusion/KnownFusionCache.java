package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.EXON_DEL_DUP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class KnownFusionCache
{
    private final List<KnownFusionData> mData;
    private final Map<KnownFusionType,List<KnownFusionData>> mDataByType;

    // cached since so commonly checked
    private final List<KnownFusionData> mKnownPairData;
    private final List<KnownFusionData> mIgRegionData;
    private final List<KnownFusionData> mHighImpactPromiscuousData;

    private boolean mHasValidData;

    public static final String KNOWN_FUSIONS_FILE = "known_fusion_file";
    public static final String KNOWN_FUSIONS_FILE_DESC = "Known fusion reference data file";
    private static final String FILE_DELIMITER = ",";

    public static final Logger KF_LOGGER = LogManager.getLogger(KnownFusionCache.class);

    public KnownFusionCache()
    {
        mData = Lists.newArrayList();
        mDataByType = Maps.newHashMap();
        mIgRegionData = Lists.newArrayList();
        mKnownPairData = Lists.newArrayList();
        mHighImpactPromiscuousData = Lists.newArrayList();
        mHasValidData = true;

        // initialise to avoid having to check for null
        Arrays.stream(KnownFusionType.values()).filter(x -> x != NONE).forEach(x -> mDataByType.put(x, Lists.newArrayList()));
    }

    public static void addKnownFusionFileOption(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(KNOWN_FUSIONS_FILE, false, KNOWN_FUSIONS_FILE_DESC);
    }

    public boolean hasValidData() { return mHasValidData; }
    public List<KnownFusionData> getData() { return mData; }
    public List<KnownFusionData> knownPairData() { return mKnownPairData; }
    public List<KnownFusionData> getDataByType(final KnownFusionType type) { return mDataByType.get(type); }

    public boolean hasKnownFusion(final String fiveGene, final String threeGene)
    {
        return mDataByType.get(KNOWN_PAIR).stream().anyMatch(x -> x.FiveGene.equals(fiveGene) && x.ThreeGene.equals(threeGene));
    }

    public boolean hasKnownUnmappable3Fusion(final String fiveGene, final String threeGene)
    {
        return mKnownPairData.stream()
                .filter(x -> !x.getThreeGeneAltRegions().isEmpty())
                .anyMatch(x -> x.FiveGene.equals(fiveGene) && x.ThreeGene.equals(threeGene));
    }

    public boolean isSingleBreakendCandidate(final String geneName, boolean isUpstream)
    {
        if(mKnownPairData.stream()
                .anyMatch(x -> (isUpstream && x.FiveGene.equals(geneName)) || (!isUpstream && x.ThreeGene.equals(geneName))))
        {
            return true;
        }

        if(mDataByType.get(IG_KNOWN_PAIR).stream()
                .filter(x -> !x.getThreeGeneAltRegions().isEmpty())
                .anyMatch(x -> !isUpstream && x.ThreeGene.equals(geneName)))
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

    public boolean hasAnyIgFusion(final String gene)
    {
        return mDataByType.get(IG_KNOWN_PAIR).stream().anyMatch(x -> x.FiveGene.equals(gene) || x.ThreeGene.equals(gene));
    }

    public boolean hasExonDelDup(final String gene)
    {
        return mDataByType.get(EXON_DEL_DUP).stream().anyMatch(x -> x.FiveGene.equals(gene) && x.ThreeGene.equals(gene));
    }

    public boolean hasKnownPairGene(final String gene)
    {
        return mKnownPairData.stream().anyMatch(x -> x.FiveGene.equals(gene) || x.ThreeGene.equals(gene));
    }

    public boolean isExonDelDupTrans(final String transName)
    {
        return mDataByType.get(EXON_DEL_DUP).stream().anyMatch(x -> x.specificExonsTransName().equals(transName));
    }

    public boolean withinPromiscuousExonRange(final KnownFusionType knownType, final String transName, int breakendExon, int fusedExon)
    {
        for(final KnownFusionData knownData : mDataByType.get(knownType))
        {
            if(!knownData.specificExonsTransName().equals(transName))
                continue;

            final int[] knownExonRange = knownType == PROMISCUOUS_5 ? knownData.fiveGeneExonRange() : knownData.threeGeneExonRange();

             if(breakendExon >= knownExonRange[SE_START] && breakendExon <= knownExonRange[SE_END]
             && fusedExon >= knownExonRange[SE_START] && fusedExon <= knownExonRange[SE_END])
             {
                 return true;
             }
        }

        return false;
    }

    public boolean withinKnownExonRanges(
            final KnownFusionType knownType, final String transName, int breakendExonUp, int fusedExonUp, int breakendExonDown, int fusedExonDown)
    {
        for(final KnownFusionData knownData : mDataByType.get(knownType))
        {
            if(!knownData.specificExonsTransName().equals(transName))
                continue;

            final int[] fiveExonRange = knownData.fiveGeneExonRange();
            final int[] threeExonRange = knownData.threeGeneExonRange();

            if(breakendExonUp >= fiveExonRange[SE_START] && breakendExonUp <= fiveExonRange[SE_END]
            && fusedExonUp >= fiveExonRange[SE_START] && fusedExonUp <= fiveExonRange[SE_END]
            && breakendExonDown >= threeExonRange[SE_START] && breakendExonDown <= threeExonRange[SE_END]
            && fusedExonDown >= threeExonRange[SE_START] && fusedExonDown <= threeExonRange[SE_END])
            {
                return true;
            }
        }

        return false;
    }

    public boolean isHighImpactPromiscuous(final KnownFusionType knownType, final String fiveGene, final String threeGene)
    {
        for(final KnownFusionData knownData : mHighImpactPromiscuousData)
        {
            if(knownData.Type != knownType)
                continue;

            if(knownData.Type == PROMISCUOUS_5 && knownData.FiveGene.equals(fiveGene))
                return true;

            if(knownData.Type == PROMISCUOUS_3 && knownData.ThreeGene.equals(threeGene))
                return true;
        }

        return false;
    }

    public boolean withinIgRegion(final String chromosome, int position)
    {
        return mIgRegionData.stream().anyMatch(x -> x.withinGeneRegion(chromosome, position));
    }

    public boolean loadFromFile(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(KNOWN_FUSIONS_FILE))
            return true;

        return loadFromFile(configBuilder.getValue(KNOWN_FUSIONS_FILE));
    }

    public boolean loadFromFile(final CommandLine cmd)
    {
        if(cmd == null || !cmd.hasOption(KNOWN_FUSIONS_FILE))
            return true;

        return loadFromFile(cmd.getOptionValue(KNOWN_FUSIONS_FILE));
    }

    public boolean loadFromFile(final String filename)
    {
        if(filename == null)
            return true;

        if(!loadFile(filename))
        {
            mHasValidData = false;
            return false;
        }

        StringJoiner refDataStr = new StringJoiner(", ");

        Arrays.stream(KnownFusionType.values())
                .filter(x -> mDataByType.containsKey(x))
                .forEach(x -> refDataStr.add(String.format("%s(%d)", x, mDataByType.get(x).size())));

        KF_LOGGER.info("loaded known fusion data: {}", refDataStr.toString());
        return true;
    }

    public void addData(final KnownFusionData data)
    {
        mData.add(data);
        mDataByType.get(data.Type).add(data);

        if(data.Type == KNOWN_PAIR)
            mKnownPairData.add(data);

        if(data.geneRegion() != null)
            mIgRegionData.add(data);

        if(data.isHighImpactPromiscuous())
            mHighImpactPromiscuousData.add(data);
    }

    public boolean loadFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            KF_LOGGER.error("file({}) not found", filename);
            return false;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
            {
                return false;
            }

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
                    KF_LOGGER.warn("file({}) invalid known fusion data will be skipped: {}", filename, data);
                }
            }
        }
        catch(IOException e)
        {
            KF_LOGGER.error("file({}) invalid known fusion data: {}", filename, e.toString());
            return false;
        }

        return true;
    }
}
