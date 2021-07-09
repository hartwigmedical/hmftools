package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_GENE_DATA_FILE;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_PROTEIN_FEATURE_DATA_FILE;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_TRANS_EXON_DATA_FILE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jooq.CSVFormat;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.tools.StringUtils;
import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class EnsemblDAO
{
    private static final String DB_URL = "ensembl_db";
    private static final String DB_USER = "ensembl_user";
    private static final String DB_PASS = "ensembl_pass";

    private static final int COORD_SYSTEM_V37 = 2;
    private static final int COORD_SYSTEM_V38 = 4;

    private DSLContext mDbContext;
    private final int mCoordSystemId;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Set<String> mGeneIds;
    private final Set<Integer> mTranscriptIds;

    public EnsemblDAO(final CommandLine cmd)
    {
        mRefGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, String.valueOf(V37)));
        mGeneIds = Sets.newHashSet();
        mTranscriptIds = Sets.newHashSet();

        mDbContext = createEnsemblDbConnection(cmd);

        if(mDbContext == null)
        {
            mCoordSystemId = -1;
            return;
        }

        mCoordSystemId = findCoordSystemId();
        GU_LOGGER.info("using coord system Id({})", mRefGenomeVersion, mCoordSystemId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - 37 (default) or 38");
        options.addOption(DB_PASS, true, "Ensembl DB password, leave out for anonymous connection");
        options.addOption(DB_URL, true, "Ensembl DB URL");
        options.addOption(DB_USER, true, "Ensembl DB username");
    }

    public boolean isValid() { return mDbContext != null && mCoordSystemId > 0; }

    private boolean ignoreTranscript(final String transName)
    {
        if(mRefGenomeVersion == V38)
            return false;

        return transName.equals("ENST00000467125"); // GOPC processed transcript which matches a ROS1 splice site
    }

    public static boolean hasDatabaseConfig(final CommandLine cmd)
    {
        return cmd.hasOption(DB_URL) && cmd.hasOption(DB_USER);
    }

    public static DSLContext createEnsemblDbConnection(final CommandLine cmd)
    {
        try
        {
            final String userName = cmd.getOptionValue(DB_USER);
            final String password = cmd.getOptionValue(DB_PASS, ""); // can be empty for an anonymous connection
            final String databaseUrl = cmd.getOptionValue(DB_URL);
            final String jdbcUrl = "jdbc:" + databaseUrl;

            System.setProperty("org.jooq.no-logo", "true");
            Connection conn = DriverManager.getConnection(jdbcUrl, userName, password);
            String catalog = conn.getCatalog();

            GU_LOGGER.info("connecting to database {}", catalog);
            return DSL.using(conn, SQLDialect.MYSQL);
        }
        catch(SQLException e)
        {
            GU_LOGGER.error("failed to connect to DB: {}", e.toString());
            return null;
        }
    }

    private int findCoordSystemId()
    {
        final String version = mRefGenomeVersion.is37() ? "GRCh37" : "GRCh38";

        final String queryStr = "select coord_system_id from coord_system c"
                + " where version = '" + version + "'"
                + " order by c.rank limit 1";

        GU_LOGGER.debug("gene query: {}", queryStr);

        Result<?> results = mDbContext.fetch(queryStr);

        for(final Record record : results)
        {
            UInteger coordSystemId = (UInteger) record.get("coord_system_id");
            return coordSystemId.intValue();
        }

        return -1;
    }

    public void writeDataCacheFiles(final String outputDir)
    {
        if(mDbContext == null)
        {
            GU_LOGGER.error("failed to establish Ensembl DB connection");
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
        GU_LOGGER.info("caching gene data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,KaryotypeBand,Synonyms");
            writer.newLine();

            String queryStr = readQueryString(Resources.getResource("sql/ensembl_gene_data.sql"));
            queryStr = queryStr.replaceAll("COORD_SYSTEM", String.valueOf(mCoordSystemId));
            // GU_LOGGER.debug("gene query: {}", queryStr);
            Result<Record> results = mDbContext.fetch(queryStr);

            final Map<String,List<EnsemblGeneData>> chrGeneMap = Maps.newHashMap();

            final Set<String> geneNames = Sets.newHashSet();

            for(final Record record : results)
            {
                String geneName = (String)record.get("GeneName");
                String geneId = (String)record.get("GeneId");

                if(geneNames.contains(geneName) || mGeneIds.contains(geneId))
                    continue;

                geneNames.add(geneName);
                mGeneIds.add(geneId);

                String chromosome = (String)record.get("Chromosome");
                UInteger geneStart = (UInteger) record.get("GeneStart");
                UInteger geneEnd = (UInteger) record.get("GeneEnd");
                Byte strand = (Byte) record.get("Strand");
                Object entrezId = record.get("EntrezId");

                if(entrezId != null)
                {
                    geneName = (String) entrezId;
                }

                EnsemblGeneData geneData = new EnsemblGeneData(
                        geneId, geneName, chromosome, strand, geneStart.intValue(), geneEnd.intValue(),
                        record.get("KaryotypeBand").toString());

                geneData.addSynonyms(record.get("Synonyms").toString());

                List<EnsemblGeneData> geneList = chrGeneMap.get(chromosome);

                if(geneList == null)
                {
                    chrGeneMap.put(chromosome, Lists.newArrayList(geneData));
                }
                else
                {
                    // add in order
                    int index = 0;
                    while(index < geneList.size())
                    {
                        if(geneData.GeneStart < geneList.get(index).GeneStart)
                            break;

                        ++index;
                    }

                    geneList.add(index, geneData);
                }
            }

            for(Map.Entry<String,List<EnsemblGeneData>> entry : chrGeneMap.entrySet())
            {
                for(EnsemblGeneData geneData : entry.getValue())
                {
                    writer.write(String.format("%s,%s,%s,%d,%d,%d,%s,%s",
                            geneData.GeneId, geneData.GeneName, geneData.Chromosome, geneData.Strand,
                            geneData.GeneStart, geneData.GeneEnd, geneData.KaryotypeBand, geneData.getSynonyms()));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
        }
    }

    private void writeTranscriptExonData(final String outputFile)
    {
        GU_LOGGER.info("caching transcript & exon data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,TransStart,TransEnd");
            writer.write(",ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd");
            writer.newLine();

            final String queryStr = readQueryString(Resources.getResource("sql/ensembl_transcript.sql"));
            // GU_LOGGER.debug("transcript query: {}", queryStr);
            Result<Record> results = mDbContext.fetch(queryStr);

            for(final Record record : results)
            {
                String transName = (String)record.get("Trans");

                if(ignoreTranscript(transName))
                    continue;

                UInteger canTransId = (UInteger) record.get("CanonicalTranscriptId");
                String geneId = (String)record.get("GeneId");
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

                if(!mGeneIds.contains(geneId))
                    continue;

                writer.write(String.format("%s,%d,%d,%d,%s,%s,%d,%d",
                        geneId, canTransId.intValue(), strand, transId.intValue(),
                        transName, record.get("BioType"), transStart.intValue(), transEnd.intValue()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%s,%s",
                        exonRank, exonStart.intValue(), exonEnd.intValue(), exonPhase, exonPhaseEnd,
                        codingStart != null ? codingStart.intValue() : "NULL",
                        codingEnd != null ? codingEnd.intValue() : "NULL"));

                writer.newLine();

                mTranscriptIds.add(transId.intValue());
            }

            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl trans-exon data file: {}", e.toString());
        }
    }

    private void writeTranscriptProteinData(final String outputFile)
    {
        GU_LOGGER.info("caching protein data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("TranscriptId,TranslationId,ProteinFeatureId,SeqStart,SeqEnd,HitDescription");
            writer.newLine();

            final String queryStr = readQueryString(Resources.getResource("sql/ensembl_protein.sql"));
            Result<Record> results = mDbContext.fetch(queryStr);

            for(final Record record : results)
            {
                UInteger transcriptId = (UInteger) record.get("TranscriptId");
                UInteger translationId = (UInteger) record.get("TranslationId");
                UInteger proteinId = (UInteger) record.get("ProteinFeatureId");
                Integer seqStart = (Integer) record.get("SeqStart");
                Integer seqEnd = (Integer) record.get("SeqEnd");

                if(!mTranscriptIds.contains(transcriptId.intValue()))
                    continue;

                writer.write(String.format("%d,%d,%d,%d,%d,%s",
                        transcriptId.intValue(), translationId.intValue(), proteinId.intValue(),
                        seqStart, seqEnd, record.get("HitDescription")));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl trans-protein data file: {}", e.toString());
        }
    }

}
