package com.hartwig.hmftools.svanalysis.gene;

import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import static org.ensembl.database.homo_sapiens_core.tables.CoordSystem.COORD_SYSTEM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class EnsemblDAO
{
    private static String DB_URL = "ensembl_db";
    private static String DB_USER = "ensembl_user";
    private static String DB_PASS = "ensembl_pass";

    public static String ENSEMBL_GENE_DATA_FILE = "ensembl_gene_data.csv";
    public static String ENSEMBL_TRANS_EXON_DATA_FILE = "ensembl_trans_exon_data.csv";
    public static String ENSEMBL_PROTEIN_FEATURE_DATA_FILE = "ensembl_protein_features.csv";

    private DSLContext mDbContext;
    private UInteger mCoordSystemId;

    private static final Logger LOGGER = LogManager.getLogger(EnsemblDAO.class);

    public EnsemblDAO(final CommandLine cmd)
    {
        mCoordSystemId = null;

        if(connectDB(cmd))
            mCoordSystemId = findCoordSystemId();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DB_PASS, true, "Ensembl DB password");
        options.addOption(DB_URL, true, "Ensembl DB URL");
        options.addOption(DB_USER, true, "Ensembl DB username");
    }

    private boolean connectDB(final CommandLine cmd)
    {
        mDbContext = null;

        try
        {
            final String userName = cmd.getOptionValue(DB_USER);
            final String password = cmd.getOptionValue(DB_PASS);
            final String databaseUrl = cmd.getOptionValue(DB_URL);
            final String jdbcUrl = "jdbc:" + databaseUrl;
            DatabaseAccess dbAccess = new DatabaseAccess(userName, password, jdbcUrl);
            mDbContext = dbAccess.context();
        }
        catch(SQLException e)
        {
            return false;
        }

        return true;
    }

    @NotNull
    private UInteger findCoordSystemId()
    {
        return mDbContext.select(COORD_SYSTEM.COORD_SYSTEM_ID)
                .from(COORD_SYSTEM)
                .where(COORD_SYSTEM.VERSION.eq("GRCh37"))
                .orderBy(COORD_SYSTEM.RANK)
                .limit(1)
                .fetchOne()
                .value1();
    }

    public void writeDataCacheFiles(final String outputDir)
    {
        if(mDbContext == null)
        {
            LOGGER.error("failed to establish Ensembl DB connection");
            return;
        }

        String outputFilename = outputDir;

        if (!outputFilename.endsWith(File.separator))
            outputFilename += File.separator;

        writeGeneData(outputFilename + ENSEMBL_GENE_DATA_FILE);
        writeTranscriptExonData(outputFilename + ENSEMBL_TRANS_EXON_DATA_FILE);
        writeTranscriptProteinData(outputFilename + ENSEMBL_PROTEIN_FEATURE_DATA_FILE);
    }

    private void writeGeneData(final String outputFile)
    {
        LOGGER.info("caching gene data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,EntrezIds,KaryotypeBand,Synonyms");
            writer.newLine();

            Result<?> results = queryAllGeneData();

            for(final Record record : results)
            {
                UInteger geneStart = (UInteger) record.get("GeneStart");
                UInteger geneEnd = (UInteger) record.get("GeneEnd");
                Byte strand = (Byte) record.get("Strand");

                writer.write(String.format("%s,%s,%s,%d,%d,%d,%s,%s,%s",
                        record.get("GeneId"), record.get("GeneName"), record.get("Chromosome"),
                        strand.intValue(), geneStart.intValue(), geneEnd.intValue(),
                        record.get("EntrezIds"), record.get("KaryotypeBand"), record.get("Synonyms")));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
        }
    }

    private Result<?> queryAllGeneData()
    {
        final String queryStr = "select gene.stable_id as GeneId,  display_xref.display_label as GeneName, seq_region.name as Chromosome,"
                + "gene.seq_region_strand as Strand, gene.seq_region_start as GeneStart, gene.seq_region_end as GeneEnd,"
                + " GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ';') as EntrezIds,"
                + " GROUP_CONCAT(DISTINCT karyotype.band ORDER BY karyotype.band SEPARATOR '-') as KaryotypeBand,"
                + " GROUP_CONCAT(DISTINCT syn_xref.dbprimary_acc ORDER BY syn_xref.dbprimary_acc SEPARATOR ';') as Synonyms"
                + " from gene"
                + " inner join object_xref as ox on gene.gene_id = ox.ensembl_id and ox.ensembl_object_type = 'GENE'"
                + " inner join xref as display_xref on display_xref.xref_id = gene.display_xref_id"
                + " inner join karyotype on gene.seq_region_id = karyotype.seq_region_id"
                + " inner join seq_region on gene.seq_region_id = seq_region.seq_region_id"
                + " left join xref as entrez_xref on (entrez_xref.xref_id = ox.xref_id and entrez_xref.external_db_id = 1300)"
                + " inner join xref as syn_xref on syn_xref.xref_id = ox.xref_id"
                + " where seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT')"
                + " and ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)"
                + " or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))"
                + " and seq_region.coord_system_id in (select coord_system_id from coord_system where version = 'GRCh37' and name = 'chromosome')"
                + " group by Chromosome, GeneStart, GeneEnd, GeneId, GeneName, Strand"
                + " order by Chromosome, GeneStart;";

        return mDbContext.fetch(queryStr);
    }

    // GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,EntrezIds,KaryotypeBand,Synonyms
    private static int GD_ID = 0;
    private static int GD_NAME = 1;
    private static int GD_CHR = 2;
    private static int GD_STRAND = 3;
    private static int GD_START = 4;
    private static int GD_END = 5;
    private static int GD_ENTREZ = 6;
    private static int GD_BAND = 7;
    private static int GD_SYN = 8;

    public static boolean loadEnsemblGeneData(final String dataPath, Map<String,
            List<EnsemblGeneData>> chrGeneDataMap, Map<String, List<EnsemblGeneData>> chrReverseGeneDataMap)
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
            List<EnsemblGeneData> reverseGeneList = null;
            String currentChr = "";
            int geneCount = 0;

            while (line != null)
            {
                String[] items = line.split(",");

                final String geneId = items[GD_ID];
                final String chromosome = items[GD_CHR];

                EnsemblGeneData geneData = new EnsemblGeneData(
                        geneId, items[GD_NAME], chromosome, Byte.parseByte(items[GD_STRAND]),
                        Long.parseLong(items[GD_START]), Long.parseLong(items[GD_END]),
                        items[GD_ENTREZ], items[GD_BAND], items[GD_SYN]);

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;
                    geneList = chrGeneDataMap.get(chromosome);

                    if(geneList == null)
                    {
                        geneList = Lists.newArrayList();
                        reverseGeneList = Lists.newArrayList();
                        chrGeneDataMap.put(chromosome, geneList);
                        chrReverseGeneDataMap.put(chromosome, reverseGeneList);
                    }
                }

                // genes are already sorted by GeneStart
                geneData.setListIndex(geneList.size());
                geneList.add(geneData);
                ++geneCount;

                // but also create a list ordered by GeneEnd
                int index = 0;
                for(; index < reverseGeneList.size(); ++index)
                {
                    final EnsemblGeneData rgd = reverseGeneList.get(index);

                    if(geneData.GeneEnd < rgd.GeneEnd)
                        break;
                }

                reverseGeneList.add(index, geneData);

                line = fileReader.readLine();
            }

            // set indicies for the reverse list
            for(Map.Entry<String, List<EnsemblGeneData>> entry : chrReverseGeneDataMap.entrySet())
            {
                final List<EnsemblGeneData> geneDataList = entry.getValue();
                for(int index = 0; index < geneDataList.size(); ++index)
                {
                    geneDataList.get(index).setReverseListIndex(index);
                }
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

    private void writeTranscriptExonData(final String outputFile)
    {
        LOGGER.info("caching transcript & exon data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,CanonicalTranscriptId,Strand,TransId,Trans,BioType,TransStart,TransEnd");
            writer.write(",ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd");
            writer.newLine();

            Result<?> results = queryAllTranscriptExonData();

            for(final Record record : results)
            {
                UInteger canTransId = (UInteger) record.get("CanonicalTranscriptId");
                UInteger transId = (UInteger) record.get("TransId");
                Byte strand = (Byte) record.get("Strand");
                UInteger transStart = (UInteger) record.get("TransStart");
                UInteger transEnd = (UInteger) record.get("TransEnd");
                Integer exonRank = (Integer) record.get("ExonRank");
                UInteger exonStart = (UInteger) record.get("ExonStart");
                UInteger exonEnd = (UInteger) record.get("ExonEnd");
                Byte exonPhase = (Byte) record.get("ExonPhase");
                Byte exonPhaseEnd = (Byte) record.get("ExonEndPhase");
                ULong codingStart = (ULong) record.get("CodingStart");
                ULong codingEnd = (ULong) record.get("CodingEnd");

                writer.write(String.format("%s,%d,%d,%d,%s,%s,%d,%d",
                        record.get("GeneId"), canTransId.intValue(), strand, transId.intValue(),
                        record.get("Trans"), record.get("BioType"), transStart.intValue(), transEnd.intValue()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%s,%s",
                        exonRank.intValue(), exonStart.intValue(), exonEnd.intValue(), exonPhase, exonPhaseEnd,
                        codingStart != null ? codingStart.intValue() : "NULL",
                        codingEnd != null ? codingEnd.intValue() : "NULL"));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing Ensembl trans-exon data file: {}", e.toString());
        }
    }

    private Result<?> queryAllTranscriptExonData()
    {
        final String queryStr = "select q1.*,"
                + " if(Strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as CodingStart,"
                + " if(Strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as CodingEnd"
                + " from ("
                + " select g.stable_id As GeneId, g.canonical_transcript_id as CanonicalTranscriptId,"
                + " t.seq_region_strand as Strand, t.transcript_id as TransId, t.stable_id as Trans, t.biotype as BioType,"
                + " t.seq_region_start as TransStart, t.seq_region_end as TransEnd,"
                + " et.rank as ExonRank, e.seq_region_start as ExonStart, e.seq_region_end as ExonEnd, e.phase as ExonPhase, e.end_phase as ExonEndPhase"
                + " from transcript as t, exon as e, exon_transcript as et, gene as g, xref as x"
                + " where t.transcript_id = et.transcript_id and e.exon_id = et.exon_id and g.display_xref_id = x.xref_id"
                + " and t.gene_id = g.gene_id"
                + " ) as q1"
                + " left join translation tl on tl.transcript_id = TransId"
                + " left join exon cs on cs.exon_id = tl.start_exon_id"
                + " left join exon ce on ce.exon_id = tl.end_exon_id"
                + " order by GeneId, TransId, ExonStart";

        return mDbContext.fetch(queryStr);
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

    public static boolean loadTranscriptExonData(final String dataPath, Map<String, List<TranscriptExonData>> geneTransExonDataMap,
            List<String> restrictedGeneIds)
    {
        String filename = dataPath;

        filename += ENSEMBL_TRANS_EXON_DATA_FILE;

        if (!Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl gene-exon data file({})", filename);
                return false;
            }

            int exonCount = 0;
            String currentGene = "";
            String lastSkippedGeneId = "";
            List<TranscriptExonData> transExonDataList = null;

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String geneId = items[TE_GENE_ID];

                if(lastSkippedGeneId.equals(geneId) || (!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneId)))
                {
                    lastSkippedGeneId = geneId;
                    line = fileReader.readLine();
                    continue;
                }

                if(!geneId.equals(currentGene))
                {
                    currentGene = geneId;
                    transExonDataList = Lists.newArrayList();
                    geneTransExonDataMap.put(geneId, transExonDataList);
                }

                // Gene,CanonicalTranscriptId,Strand,TransId,Trans,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,
                // ExonPhase,ExonEndPhase,CodingStart,CodingEnd

                Long codingStart = !items[TE_CODING_START].equalsIgnoreCase("NULL") ? Long.parseLong(items[TE_CODING_START]) : null;
                Long codingEnd = !items[TE_CODING_END].equalsIgnoreCase("NULL") ? Long.parseLong(items[TE_CODING_END]) : null;
                int transId = Integer.parseInt(items[TE_TRANS_ID]);
                int canonicalTransId = Integer.parseInt(items[TE_CANONICAL]);

                TranscriptExonData exonData = new TranscriptExonData(
                        geneId, items[TE_TRANS_NAME], transId, transId == canonicalTransId, Byte.parseByte(items[TE_STRAND]),
                        Long.parseLong(items[TE_TRANS_START]), Long.parseLong(items[TE_TRANS_END]),
                        Long.parseLong(items[TE_EXON_START]), Long.parseLong(items[TE_EXON_END]),
                        Integer.parseInt(items[TE_EXON_RANK]), Integer.parseInt(items[TE_PHASE]), Integer.parseInt(items[TE_PHASE_END]),
                        codingStart, codingEnd, items[TE_BIOTYPE]);

                transExonDataList.add(exonData);
                ++exonCount;

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} gene records, {} exons", geneTransExonDataMap.size(), exonCount);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load gene transcript exon data({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private void writeTranscriptProteinData(final String outputFile)
    {
        LOGGER.info("caching protein data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("TranscriptId,TranslationId,ProteinFeatureId,SeqStart,SeqEnd,HitDescription");
            writer.newLine();

            Result<?> results = queryAllProteinData();

            for(final Record record : results)
            {
                UInteger transcriptId = (UInteger) record.get("TranscriptId");
                UInteger translationId = (UInteger) record.get("TranslationId");
                UInteger proteinId = (UInteger) record.get("ProteinFeatureId");
                Integer seqStart = (Integer) record.get("SeqStart");
                Integer seqEnd = (Integer) record.get("SeqEnd");

                writer.write(String.format("%d,%d,%d,%d,%d,%s",
                        transcriptId.intValue(), translationId.intValue(), proteinId.intValue(),
                        seqStart.intValue(), seqEnd.intValue(), record.get("HitDescription")));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing Ensembl trans-protein data file: {}", e.toString());
        }
    }

    private Result<?> queryAllProteinData()
    {
        final String queryStr = "select tl.transcript_id as TranscriptId, tl.translation_id as TranslationId, protein_feature_id as ProteinFeatureId,"
                + " pf.seq_start as SeqStart, pf.seq_end as SeqEnd, hit_description as HitDescription"
                + " from protein_feature pf, analysis_description ad, translation tl, transcript t"
                + " where pf.analysis_id = ad.analysis_id and pf.translation_id = tl.translation_id and t.transcript_id = tl.transcript_id"
                + " and display_label = 'PROSITE profiles'"
                + " order by tl.transcript_id, tl.translation_id, pf.seq_start";

        return mDbContext.fetch(queryStr);
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



}
