package com.hartwig.hmftools.geneutils.fusion;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.OVERRIDE_DOWN_DISTANCE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.OVERRIDE_IG_RANGE;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR_DESC;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.createOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.getEnsemblDirectory;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class GenerateFusionFiles
{
    private final String mKnownFusionDbFile;
    private final String mResourceRepoDir;
    private final String mOutputDir;

    private static final String KNOWN_FUSION_DB_FILE = "known_fusion_db_file";
    private static final int PRE_GENE_BUFFER = 10000;

    public GenerateFusionFiles(final ConfigBuilder configBuilder)
    {
        GU_LOGGER.info("starting known fusion file generation");

        mKnownFusionDbFile = configBuilder.getValue(KNOWN_FUSION_DB_FILE);
        mResourceRepoDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_REPO_DIR));
        mOutputDir = parseOutputDir(configBuilder);
    }

    public void run()
    {
        // step 1: load fusion knowledge database file
        List<FusionRefData> fusionRefData = loadFusionRefData();

        if(fusionRefData.isEmpty())
            System.exit(1);

        // step 2: check for duplicates
        for(int i = 0; i < fusionRefData.size() - 1; ++i)
        {
            FusionRefData ref1 = fusionRefData.get(i);

            for(int j = i + 1; j < fusionRefData.size(); ++j)
            {
                FusionRefData ref2 = fusionRefData.get(j);

                if(ref1.isDuplicate(ref2))
                {
                    GU_LOGGER.error("duplicate fusion entry: {}", ref1.toString());
                    System.exit(1);
                }
            }
        }

        createOutputDir(mOutputDir);
        createFusionFiles(RefGenomeVersion.V37, fusionRefData);
        createFusionFiles(RefGenomeVersion.V38, fusionRefData);

        GU_LOGGER.info("fusion reference file generation complete");
    }

    private void createFusionFiles(final RefGenomeVersion refGenomeVersion, final List<FusionRefData> fusionRefData)
    {
        // step 3: load Ensembl data cache files
        String ensemblDir = getEnsemblDirectory(refGenomeVersion, mResourceRepoDir);

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);
        ensemblDataCache.createGeneNameIdMap();

        // step 4: check all genes are in Ensembl
        for(FusionRefData fusion : fusionRefData)
        {
            if(!fusion.FiveGene.isEmpty() && !fusion.FiveGene.startsWith("IG"))
            {
                if(ensemblDataCache.getGeneDataByName(fusion.FiveGene) == null)
                {
                    GU_LOGGER.error("fusion({}) five-prime gene not present in Ensembl", fusion);
                    System.exit(1);
                }

            }

            if(!fusion.ThreeGene.isEmpty() && !fusion.ThreeGene.equals("C19MC"))
            {
                if(ensemblDataCache.getGeneDataByName(fusion.ThreeGene) == null)
                {
                    GU_LOGGER.error("fusion({}) three-prime gene not present in Ensembl", fusion);
                    System.exit(1);
                }
            }
        }

        // step 5: write known_fusion_data for each version
        writeKnownFusionFiles(refGenomeVersion, fusionRefData);

        // step 6: write known fusion BED files (eg for Gripss and SvPrep)
        writeFusionBedFiles(refGenomeVersion, fusionRefData, ensemblDataCache);
    }

    private void writeKnownFusionFiles(final RefGenomeVersion refGenomeVersion, final List<FusionRefData> fusionRefData)
    {
        try
        {
            String fusionFilename = refGenomeVersion.addVersionToFilePath(mOutputDir + "known_fusion_data.csv");

            BufferedWriter writer = createBufferedWriter(fusionFilename, false);

            StringJoiner header = new StringJoiner(KnownFusionData.FILE_DELIM);
            header.add(KnownFusionData.FLD_TYPE);
            header.add(KnownFusionData.FLD_FIVE_GENE);
            header.add(KnownFusionData.FLD_THREE_GENE);
            header.add(KnownFusionData.FLD_CANCER_TYPES);
            header.add(KnownFusionData.FLD_PUB_MED);
            header.add(KnownFusionData.FLD_KNOWN_EXON_TRANS);
            header.add(KnownFusionData.FLD_KNOWN_EXON_UP_RANGE);
            header.add(KnownFusionData.FLD_KNOWN_EXON_DOWN_RANGE);
            header.add(KnownFusionData.FLD_HIGH_IMPACT_PROM);
            header.add(KnownFusionData.FLD_OVERRIDES);

            writer.write(header.toString());
            writer.newLine();

            for(FusionRefData fusion : fusionRefData)
            {
                StringJoiner fusionData = new StringJoiner(KnownFusionData.FILE_DELIM);
                fusionData.add(fusion.Type.toString());
                fusionData.add(fusion.FiveGene);
                fusionData.add(fusion.ThreeGene);
                fusionData.add(fusion.CancerTypes);
                fusionData.add(fusion.PubMedId);

                if(refGenomeVersion.is37())
                {
                    fusionData.add(fusion.KnownExonTranscript);
                    fusionData.add(fusion.KnownExonUpRange);
                    fusionData.add(fusion.KnownExonDownRange);
                }
                else
                {
                    fusionData.add(fusion.KnownExonTranscriptRef38);
                    fusionData.add(fusion.KnownExonUpRangeRef38);
                    fusionData.add(fusion.KnownExonDownRangeRef38);
                }

                fusionData.add(fusion.HighImpactPromiscuous);

                if(refGenomeVersion.is37())
                    fusionData.add(fusion.Overrides);
                else
                    fusionData.add(fusion.OverridesRef38);

                writer.write(fusionData.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write known fusion data file: {}", e.toString());
        }
    }

    private static class FusionBedData implements Comparable<FusionBedData>
    {
        public final String Name;
        public final String ChrUp;
        public final String ChrDown;
        public final byte StrandUp;
        public final byte StrandDown;
        public final int PositionUpStart;
        public final int PositionUpEnd;
        public final int PositionDownStart;
        public final int PositionDownEnd;

        private final boolean mUpIsStart;

        public FusionBedData(
                final String name, final String chrUp, final String chrDown, final byte strandUp, final byte strandDown,
                final int positionUpStart, final int positionUpEnd, final int positionDownStart, final int positionDownEnd)
        {
            Name = name;
            ChrUp = chrUp;
            ChrDown = chrDown;
            StrandUp = strandUp;
            StrandDown = strandDown;
            PositionUpStart = positionUpStart;
            PositionUpEnd = positionUpEnd;
            PositionDownStart = positionDownStart;
            PositionDownEnd = positionDownEnd;

            if(lowerChromosome(ChrUp, ChrDown))
                mUpIsStart = true;
            else if(ChrUp.equals(ChrDown))
                mUpIsStart = adjustedUpStart() < adjustedDownStart();
            else
                mUpIsStart = false;
        }

        // methods for sorting
        public int adjustedUpStart() { return StrandUp == POS_STRAND ? PositionUpStart - PRE_GENE_BUFFER - 1 : PositionUpStart - 1; }
        public int adjustedUpEnd() { return StrandUp == POS_STRAND ? PositionUpEnd : PositionUpEnd + PRE_GENE_BUFFER; }
        public int adjustedDownStart() { return StrandDown == POS_STRAND ? PositionDownStart - PRE_GENE_BUFFER - 1 : PositionDownStart - 1; }
        public int adjustedDownEnd() { return StrandDown == POS_STRAND ? PositionDownEnd : PositionDownEnd + PRE_GENE_BUFFER; }

        public char strandUpChar() { return StrandUp == POS_STRAND ? '+' : '-'; }

        // REVERSE STRAND2 since for the downstream genes the orientation is opposite to upstream (ie +ve strand = -ve orientation and vice versa)
        public char strandDownChar() { return StrandDown == POS_STRAND ? '-' : '+'; }

        // ordered for the BED
        public String chrStart() { return mUpIsStart ? ChrUp : ChrDown; }
        public String chrEnd() { return !mUpIsStart ? ChrUp : ChrDown; }
        public int posStartStart() { return mUpIsStart ? adjustedUpStart() : adjustedDownStart(); }
        public int posStartEnd() { return mUpIsStart ? adjustedUpEnd() : adjustedDownEnd(); }
        public int posEndStart() { return !mUpIsStart ? adjustedUpStart() : adjustedDownStart(); }
        public int posEndEnd() { return !mUpIsStart ? adjustedUpEnd() : adjustedDownEnd(); }
        public char strandStart() { return mUpIsStart ? strandUpChar() : strandDownChar(); }
        public char strandEnd() { return !mUpIsStart ? strandUpChar() : strandDownChar(); }

        public String toString()
        {
            return format("%s start(%s:%d-%d) end(%s:%d-%d) upIsStart(%s)",
                    Name, chrStart(), posStartStart(), posStartEnd(), chrEnd(), posEndStart(), posEndEnd(), mUpIsStart);
        }

        @Override
        public int compareTo(final FusionBedData other)
        {
            if(lowerChromosome(chrStart(), other.chrStart()))
            {
                return -1;
            }
            else if(chrStart().equals(other.chrStart()))
            {
                if(posStartStart() == other.posStartStart())
                    return Name.compareTo(other.Name);

                return posStartStart() < other.posStartStart() ? -1 : 1;
            }
            else
            {
                return 1;
            }
        }
    }

    private void writeFusionBedFiles(
            final RefGenomeVersion refGenomeVersion, final List<FusionRefData> fusionRefData, final EnsemblDataCache ensemblDataCache)
    {
        List<FusionBedData> fusionBedDataList = Lists.newArrayList();

        for(FusionRefData fusion : fusionRefData)
        {
            addBedEntries(refGenomeVersion, ensemblDataCache, fusion, fusionBedDataList);
        }

        // sort and then write
        Collections.sort(fusionBedDataList);

        try
        {
            String bedFilename = refGenomeVersion.addVersionToFilePath(mOutputDir + "known_fusions.bedpe");

            BufferedWriter writer = createBufferedWriter(bedFilename, false);

            for(FusionBedData fusionBedData : fusionBedDataList)
            {
                StringJoiner fusionData = new StringJoiner(TSV_DELIM);

                fusionData.add(refGenomeVersion.versionedChromosome(fusionBedData.chrStart()));
                fusionData.add(String.valueOf(fusionBedData.posStartStart()));
                fusionData.add(String.valueOf(fusionBedData.posStartEnd()));
                fusionData.add(refGenomeVersion.versionedChromosome(fusionBedData.chrEnd()));
                fusionData.add(String.valueOf(fusionBedData.posEndStart()));
                fusionData.add(String.valueOf(fusionBedData.posEndEnd()));
                fusionData.add(fusionBedData.Name);
                fusionData.add("0"); // score isn't used
                fusionData.add(String.valueOf(fusionBedData.strandStart()));
                fusionData.add(String.valueOf(fusionBedData.strandEnd()));

                writer.write(fusionData.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write known fusion data file: {}", e.toString());
        }
    }

    private void addBedEntries(
            final RefGenomeVersion refGenomeVersion, final EnsemblDataCache ensemblDataCache, final FusionRefData fusion,
            final List<FusionBedData> bedEntries)
    {
        if(fusion.Type != KNOWN_PAIR && fusion.Type != IG_KNOWN_PAIR)
            return;

        String chrUp;
        String chrDown;
        int posUpStart;
        int posUpEnd;
        int posDownStart;
        int posDownEnd;
        byte strandUp;
        byte strandDown;

        String name = format("%s-%s", fusion.FiveGene, fusion.ThreeGene);

        GeneData threeGeneData = ensemblDataCache.getGeneDataByName(fusion.ThreeGene);

        chrDown = threeGeneData.Chromosome;
        posDownStart = threeGeneData.GeneStart;
        posDownEnd = threeGeneData.GeneEnd;
        strandDown = threeGeneData.Strand;

        if(fusion.Type == KNOWN_PAIR)
        {
            GeneData fiveGeneData = ensemblDataCache.getGeneDataByName(fusion.FiveGene);

            chrUp = fiveGeneData.Chromosome;
            posUpStart = fiveGeneData.GeneStart;
            posUpEnd = fiveGeneData.GeneEnd;
            strandUp = fiveGeneData.Strand;

            bedEntries.add(new FusionBedData(
                    name, chrUp, chrDown, strandUp, strandDown, posUpStart, posUpEnd, posDownStart, posDownEnd));
        }
        else
        {
            String overrides = refGenomeVersion.is37() ? fusion.Overrides : fusion.OverridesRef38;
            String[] overridesStr = overrides.split(" ", -1);

            String igRangeStr = overridesStr[0].replaceAll(OVERRIDE_IG_RANGE + "=", "");
            String[] igRangeItems = igRangeStr.split(ITEM_DELIM);
            strandUp = Byte.parseByte(igRangeItems[0]);
            byte upStrandReversed = strandUp == POS_STRAND ? NEG_STRAND : POS_STRAND;
            chrUp = igRangeItems[1];
            posUpStart = Integer.parseInt(igRangeItems[2]);
            posUpEnd = Integer.parseInt(igRangeItems[3]);

            bedEntries.add(new FusionBedData(
                    name, chrUp, chrDown, strandUp, strandDown, posUpStart, posUpEnd, posDownStart, posDownEnd));

            // add a second entry with the IG gene reversed
            bedEntries.add(new FusionBedData(
                    name, chrUp, chrDown, upStrandReversed, strandDown, posUpStart, posUpEnd, posDownStart, posDownEnd));

            // add another entry for the downstream location if present
            if(overridesStr.length > 1 && overridesStr[1].contains(OVERRIDE_DOWN_DISTANCE))
            {
                String downstreamStr = overridesStr[1];
                downstreamStr = downstreamStr.replaceAll(OVERRIDE_DOWN_DISTANCE + "=", "");
                int downstreamDistance = Integer.parseInt(downstreamStr);

                posDownStart = strandDown == POS_STRAND ?
                        threeGeneData.GeneEnd : threeGeneData.GeneStart - downstreamDistance;

                posDownEnd = strandDown == POS_STRAND ?
                        threeGeneData.GeneEnd + downstreamDistance : threeGeneData.GeneStart;

                strandDown = strandDown == POS_STRAND ? NEG_STRAND : POS_STRAND;

                bedEntries.add(new FusionBedData(
                        name, chrUp, chrDown, strandUp, strandDown, posUpStart, posUpEnd, posDownStart, posDownEnd));

                bedEntries.add(new FusionBedData(
                        name, chrUp, chrDown, upStrandReversed, strandDown, posUpStart, posUpEnd, posDownStart, posDownEnd));
            }
        }
    }

    private List<FusionRefData> loadFusionRefData()
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(mKnownFusionDbFile));
            String header = lines.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            List<FusionRefData> fusions = Lists.newArrayList();

            for(int i = 1; i < lines.size(); ++i)
            {
                String line = lines.get(i);
                String[] values = line.split(TSV_DELIM, -1);

                if(values.length < 13)
                {
                    GU_LOGGER.error("entry({}) missing data: {}", i, line);
                    return Collections.emptyList();
                }

                KnownFusionType type = KnownFusionType.valueOf(values[fieldsIndexMap.get(KnownFusionData.FLD_TYPE)]);

                fusions.add(new FusionRefData(
                        type,
                        values[fieldsIndexMap.get(KnownFusionData.FLD_FIVE_GENE)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_THREE_GENE)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_CANCER_TYPES)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_PUB_MED)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_KNOWN_EXON_TRANS)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_KNOWN_EXON_UP_RANGE)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_KNOWN_EXON_DOWN_RANGE)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_HIGH_IMPACT_PROM)],
                        values[fieldsIndexMap.get(KnownFusionData.FLD_OVERRIDES)],
                        values[fieldsIndexMap.get("KnownExonTranscriptRef38")],
                        values[fieldsIndexMap.get("KnownExonUpRangeRef38")],
                        values[fieldsIndexMap.get("KnownExonDownRangeRef38")],
                        values[fieldsIndexMap.get("OverridesRef38")]));
            }

            GU_LOGGER.info("loaded {} fusion entries from file: {}", fusions.size(), mKnownFusionDbFile);
            return fusions;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to read fusion file: {}", e.toString());
            return Collections.emptyList();
        }
    }

    private static class FusionRefData
    {
        public final KnownFusionType Type;
        public final String FiveGene;
        public final String ThreeGene;
        public final String CancerTypes;
        public final String PubMedId;
        public final String KnownExonTranscript;
        public final String KnownExonUpRange;
        public final String KnownExonDownRange;
        public final String HighImpactPromiscuous;
        public final String Overrides;
        public final String KnownExonTranscriptRef38;
        public final String KnownExonUpRangeRef38;
        public final String KnownExonDownRangeRef38;
        public final String OverridesRef38;

        public FusionRefData(
                final KnownFusionType type, final String fiveGene, final String threeGene, final String cancerTypes, final String pubMedId,
                final String knownExonTranscript, final String knownExonUpRange, final String knownExonDownRange,
                final String highImpactPromiscuous, final String overrides, final String knownExonTranscriptRef38,
                final String knownExonUpRangeRef38, final String knownExonDownRangeRef38, final String overridesRef38)
        {
            Type = type;
            FiveGene = fiveGene;
            ThreeGene = threeGene;
            CancerTypes = cancerTypes;
            PubMedId = pubMedId;
            KnownExonTranscript = knownExonTranscript;
            KnownExonUpRange = knownExonUpRange;
            KnownExonDownRange = knownExonDownRange;
            HighImpactPromiscuous = highImpactPromiscuous;
            Overrides = overrides;
            KnownExonTranscriptRef38 = knownExonTranscriptRef38;
            KnownExonUpRangeRef38 = knownExonUpRangeRef38;
            KnownExonDownRangeRef38 = knownExonDownRangeRef38;
            OverridesRef38 = overridesRef38;
        }

        public boolean isDuplicate(final FusionRefData other)
        {
            if(Type != other.Type)
                return false;

            if(!FiveGene.equals(other.FiveGene) || !ThreeGene.equals(other.ThreeGene))
                return false;

            boolean equalKnownExonRange37 =
                    KnownExonTranscript.equals(other.KnownExonTranscript) && KnownExonUpRange.equals(other.KnownExonUpRange)
                            && KnownExonDownRange.equals(other.KnownExonDownRange);
            boolean equalKnownExonRange38 = KnownExonTranscriptRef38.equals(other.KnownExonTranscriptRef38)
                    && KnownExonUpRangeRef38.equals(other.KnownExonUpRangeRef38)
                    && KnownExonDownRangeRef38.equals(other.KnownExonDownRangeRef38);
            return equalKnownExonRange37 || equalKnownExonRange38;
        }

        public String toString()
        {
            return format("%s: %s-%s", Type, FiveGene, ThreeGene);
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(KNOWN_FUSION_DB_FILE, true, "File containing the driver gene panel for 37");
        configBuilder.addPath(RESOURCE_REPO_DIR, true, RESOURCE_REPO_DIR_DESC);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateFusionFiles generator = new GenerateFusionFiles(configBuilder);
        generator.run();
    }
}
