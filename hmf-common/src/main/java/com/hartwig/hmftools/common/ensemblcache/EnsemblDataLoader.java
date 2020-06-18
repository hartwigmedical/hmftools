package com.hartwig.hmftools.common.ensemblcache;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.refGenomeChromosome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class EnsemblDataLoader
{
    public static final String ENSEMBL_GENE_DATA_FILE = "ensembl_gene_data.csv";
    public static final String ENSEMBL_TRANS_EXON_DATA_FILE = "ensembl_trans_exon_data.csv";
    public static final String ENSEMBL_TRANS_SPLICE_DATA_FILE = "ensembl_trans_splice_data.csv";
    public static final String ENSEMBL_PROTEIN_FEATURE_DATA_FILE = "ensembl_protein_features.csv";

    private static final Logger LOGGER = LogManager.getLogger(EnsemblDataLoader.class);

    // GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,EntrezIds,KaryotypeBand,Synonyms
    private static int GD_ID = 0;
    private static int GD_NAME = 1;
    private static int GD_CHR = 2;
    private static int GD_STRAND = 3;
    private static int GD_START = 4;
    private static int GD_END = 5;
    private static int GD_ENTREZ = 6; // currently unused
    private static int GD_BAND = 7;
    private static int GD_SYN = 8; // currently unused

    public static boolean loadEnsemblGeneData(final String dataPath, final List<String> restrictedGeneIds,
            final Map<String, List<EnsemblGeneData>> chrGeneDataMap, RefGenomeVersion version)
    {
        String filename = dataPath;

        filename += ENSEMBL_GENE_DATA_FILE;

        if (!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl gene data file({})", filename);
                return false;
            }

            line = fileReader.readLine(); // skip header

            List<EnsemblGeneData> geneList = null;
            String currentChr = "";
            int geneCount = 0;

            while (line != null)
            {
                String[] items = line.split(",");

                final String geneId = items[GD_ID];

                if(!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneId))
                {
                    line = fileReader.readLine();
                    continue;
                }

                final String chromosome = refGenomeChromosome(items[GD_CHR], version);

                EnsemblGeneData geneData = new EnsemblGeneData(
                        geneId, items[GD_NAME], chromosome, Byte.parseByte(items[GD_STRAND]),
                        Integer.parseInt(items[GD_START]), Integer.parseInt(items[GD_END]), items[GD_BAND]);

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;
                    geneList = chrGeneDataMap.get(chromosome);

                    if(geneList == null)
                    {
                        geneList = Lists.newArrayList();
                        chrGeneDataMap.put(chromosome, geneList);
                    }
                }

                geneList.add(geneData);
                ++geneCount;

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} gene records", geneCount);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load Ensembl gene ({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }


    // Gene,CanonicalTranscriptId,Strand,TransId,Trans,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,
    // ExonPhase,ExonEndPhase,CodingStart,CodingEnd
    private static int TE_GENE_ID = 0;
    private static int TE_CANONICAL = 1;
    private static int TE_STRAND = 2;
    private static int TE_TRANS_ID = 3;
    private static int TE_TRANS_NAME = 4;
    private static int TE_BIOTYPE = 5;
    private static int TE_TRANS_START = 6;
    private static int TE_TRANS_END = 7;
    private static int TE_EXON_RANK = 8;
    private static int TE_EXON_START = 9;
    private static int TE_EXON_END = 10;
    private static int TE_PHASE = 11;
    private static int TE_PHASE_END = 12;
    private static int TE_CODING_START = 13;
    private static int TE_CODING_END = 14;

    public static boolean loadTranscriptData(final String dataPath, Map<String, List<TranscriptData>> transcriptDataMap,
            final List<String> restrictedGeneIds, boolean cacheExons, boolean canonicalOnly)
    {
        String filename = dataPath;

        filename += ENSEMBL_TRANS_EXON_DATA_FILE;

        if (!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            if (line == null)
            {
                LOGGER.error("empty Ensembl gene-exon data file({})", filename);
                return false;
            }

            int exonCount = 0;
            int transcriptCount = 0;
            String currentGene = "";
            TranscriptData currentTrans = null;
            String lastSkippedGeneId = "";
            List<TranscriptData> transDataList = null;
            List<ExonData> exonDataList = null;

            while (line != null)
            {
                line = fileReader.readLine();

                if(line == null)
                    break;

                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String geneId = items[TE_GENE_ID];
                int transId = Integer.parseInt(items[TE_TRANS_ID]);

                if(lastSkippedGeneId.equals(geneId) || (!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneId)))
                {
                    lastSkippedGeneId = geneId;
                    continue;
                }

                if(!geneId.equals(currentGene))
                {
                    currentGene = geneId;
                    transDataList = Lists.newArrayList();
                    transcriptDataMap.put(geneId, transDataList);
                }

                // Gene,CanonicalTranscriptId,Strand,TransId,Trans,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,
                // ExonPhase,ExonEndPhase,CodingStart,CodingEnd

                if(currentTrans == null || currentTrans.TransId != transId)
                {
                    int canonicalTransId = Integer.parseInt(items[TE_CANONICAL]);
                    boolean isCanonical = (canonicalTransId == transId);

                    if(!isCanonical && canonicalOnly)
                        continue;

                    exonDataList = Lists.newArrayList();

                    Integer codingStart = !items[TE_CODING_START].equalsIgnoreCase("NULL") ? Integer.parseInt(items[TE_CODING_START]) : null;
                    Integer codingEnd = !items[TE_CODING_END].equalsIgnoreCase("NULL") ? Integer.parseInt(items[TE_CODING_END]) : null;

                    currentTrans = new TranscriptData(
                            transId, items[TE_TRANS_NAME], geneId, isCanonical, Byte.parseByte(items[TE_STRAND]),
                            Integer.parseInt(items[TE_TRANS_START]), Integer.parseInt(items[TE_TRANS_END]),
                            codingStart, codingEnd, items[TE_BIOTYPE]);

                    ++transcriptCount;

                    currentTrans.setExons(exonDataList);
                    transDataList.add(currentTrans);
                }

                if(cacheExons || currentTrans.IsCanonical)
                {
                    ExonData exonData = new ExonData(
                            transId, Integer.parseInt(items[TE_EXON_START]), Integer.parseInt(items[TE_EXON_END]),
                            Integer.parseInt(items[TE_EXON_RANK]), Integer.parseInt(items[TE_PHASE]), Integer.parseInt(items[TE_PHASE_END]));

                    exonDataList.add(exonData);
                    ++exonCount;
                }
            }

            LOGGER.debug("loaded {} genes with {} transcripts records and {} exons",
                    transcriptDataMap.size(), transcriptCount, exonCount);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load gene transcript exon data({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    // TranscriptId,TranslationId,ProteinFeatureId,SeqStart,SeqEnd,HitDescription
    private static int PF_TRANS_ID = 0;
    private static int PF_TRANL_ID = 1;
    private static int PF_PF_ID = 2;
    private static int PF_START = 3;
    private static int PF_END = 4;
    private static int PF_DESC = 5;

    public static boolean loadTranscriptProteinData(final String dataPath, Map<Integer, List<TranscriptProteinData>> proteinDataMap,
            List<Integer> restrictedTransIds)
    {
        String filename = dataPath;

        filename += ENSEMBL_PROTEIN_FEATURE_DATA_FILE;

        if (!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl protein feature data file({})", filename);
                return false;
            }

            int proteinCount = 0;
            int currentTransId = -1;
            List<TranscriptProteinData> transProteinDataList = null;

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                int transId = Integer.parseInt(items[PF_TRANS_ID]);

                if(!restrictedTransIds.isEmpty() && !restrictedTransIds.contains(transId))
                {
                    line = fileReader.readLine();
                    continue;
                }

                if(transId != currentTransId)
                {
                    currentTransId = transId;
                    transProteinDataList = Lists.newArrayList();
                    proteinDataMap.put(transId, transProteinDataList);
                }

                TranscriptProteinData proteinData = new TranscriptProteinData(
                        transId, Integer.parseInt(items[PF_TRANL_ID]), Integer.parseInt(items[PF_PF_ID]),
                        Integer.parseInt(items[PF_START]), Integer.parseInt(items[PF_END]), items[PF_DESC]);

                transProteinDataList.add(proteinData);
                ++proteinCount;

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} protein trans records with {} locations", proteinDataMap.size(), proteinCount);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load transcript protein features({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    // GeneId,TransId,TransName,TransStartPos,PreSpliceAcceptorPosition,Distance
    private static int TA_GENE_ID = 0;
    private static int TA_TRANS_ID = 1;
    private static int TA_TRANS_NAME = 2;
    private static int TA_TRANS_START = 3;
    private static int TA_SA_POSITION = 4;
    private static int TA_DISTANCE = 5;

    public static boolean loadTranscriptSpliceAcceptorData(
            final String dataPath, Map<Integer,Integer> transSaPositionDataMap, final List<Integer> restrictedTransIds)
    {
        String filename = dataPath;

        filename += ENSEMBL_TRANS_SPLICE_DATA_FILE;

        if (!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl trans splice acceptor data file({})", filename);
                return false;
            }

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final int transId = Integer.parseInt(items[TA_TRANS_ID]);

                if(!restrictedTransIds.isEmpty() && !restrictedTransIds.contains(transId))
                {
                    line = fileReader.readLine();
                    continue;
                }

                int saPosition = Integer.parseInt(items[TA_SA_POSITION]);

                transSaPositionDataMap.put(transId, saPosition);

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} trans splice-acceptor position records", transSaPositionDataMap.size());
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load transcript splice data({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

}
