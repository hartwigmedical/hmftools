package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class EnsemblDAO
{
    private static final String DB_URL = "ensembl_db";
    private static final String DB_USER = "ensembl_user";
    private static final String DB_PASS = "ensembl_pass";

    public static final String ENSEMBL_GENE_DATA_FILE = "ensembl_gene_data.csv";
    public static final String ENSEMBL_TRANS_EXON_DATA_FILE = "ensembl_trans_exon_data.csv";
    public static final String ENSEMBL_TRANS_SPLICE_DATA_FILE = "ensembl_trans_splice_data.csv";
    public static final String ENSEMBL_PROTEIN_FEATURE_DATA_FILE = "ensembl_protein_features.csv";

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

        if(!connectDB(cmd))
        {
            mCoordSystemId = -1;
            return;
        }

        mCoordSystemId = findCoordSystemId();
        LNX_LOGGER.info("Ref genome version({}), coord system Id({})", mRefGenomeVersion, mCoordSystemId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - 37 (default) or 38");
        options.addOption(DB_PASS, true, "Ensembl DB password");
        options.addOption(DB_URL, true, "Ensembl DB URL");
        options.addOption(DB_USER, true, "Ensembl DB username");
    }

    public boolean isValid() { return mDbContext != null && mCoordSystemId > 0; }
    public RefGenomeVersion refGenomeVersion() { return mRefGenomeVersion; }

    private boolean ignoreTranscript(final String transName)
    {
        if(mRefGenomeVersion == V38)
            return false;

        return transName.equals("ENST00000467125"); // GOPC processed transcript which matches a ROS1 splice site
    }

    public static boolean hasDatabaseConfig(final CommandLine cmd)
    {
        return cmd.hasOption(DB_URL) && cmd.hasOption(DB_USER) && cmd.hasOption(DB_PASS);
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
            LNX_LOGGER.error("failed to connect to DB: {}", e.toString());
            return false;
        }

        return true;
    }

    private int findCoordSystemId()
    {
        final String version = mRefGenomeVersion.is37() ? "GRCh37" : "GRCh38";

        final String queryStr = "select coord_system_id from coord_system c"
                + " where version = '" + version + "'"
                + " order by c.rank limit 1";

        LNX_LOGGER.debug("gene query: {}", queryStr);

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
            LNX_LOGGER.error("failed to establish Ensembl DB connection");
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
        LNX_LOGGER.info("caching gene data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,KaryotypeBand,Synonyms");
            writer.newLine();

            Result<?> results = queryAllGeneData();

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
            LNX_LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
        }
    }

    private Result<?> queryAllGeneData()
    {
        final String queryStr = "select gene.stable_id as GeneId, display_xref.display_label as GeneName, entrez_xref.display_label as EntrezId,"
                + " case when entrez_xref.external_db_id = 1100 then 1100 else 0 end as ExtDbId, seq_region.name as Chromosome,"
                + " gene.seq_region_strand as Strand, gene.seq_region_start as GeneStart, gene.seq_region_end as GeneEnd,"
                + " GROUP_CONCAT(DISTINCT entrez_xref.dbprimary_acc ORDER BY entrez_xref.dbprimary_acc SEPARATOR ';') as EntrezIds,"
                + " GROUP_CONCAT(DISTINCT karyotype.band ORDER BY karyotype.band SEPARATOR '-') as KaryotypeBand,"
                + " GROUP_CONCAT(DISTINCT syn_xref.dbprimary_acc ORDER BY syn_xref.dbprimary_acc SEPARATOR ';') as Synonyms"
                + " from gene"
                + " inner join object_xref as ox on gene.gene_id = ox.ensembl_id and ox.ensembl_object_type = 'GENE'"
                + " inner join xref as display_xref on display_xref.xref_id = gene.display_xref_id"
                + " inner join karyotype on gene.seq_region_id = karyotype.seq_region_id"
                + " inner join seq_region on gene.seq_region_id = seq_region.seq_region_id"
                + " left join xref as entrez_xref on (entrez_xref.xref_id = ox.xref_id and entrez_xref.external_db_id = 1100)"
                + " inner join xref as syn_xref on syn_xref.xref_id = ox.xref_id"
                + " where seq_region.name in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')"
                + " and ((gene.seq_region_start >= karyotype.seq_region_start and gene.seq_region_start <= karyotype.seq_region_end)"
                + " or (gene.seq_region_end >= karyotype.seq_region_start and gene.seq_region_end <= karyotype.seq_region_end))"
                + " and seq_region.coord_system_id = " + mCoordSystemId
                + " group by GeneId, GeneName, EntrezId, ExtDbId, Chromosome, Strand, GeneStart, GeneEnd"
                + " order by GeneName, ExtDbId desc;";

        LNX_LOGGER.debug("gene query: {}", queryStr);

        return mDbContext.fetch(queryStr);
    }

    private void writeTranscriptExonData(final String outputFile)
    {
        LNX_LOGGER.info("caching transcript & exon data to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,TransStart,TransEnd");
            writer.write(",ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd");
            writer.newLine();

            Result<?> results = queryAllTranscriptExonData();

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
            LNX_LOGGER.error("error writing Ensembl trans-exon data file: {}", e.toString());
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

        LNX_LOGGER.debug("transcript query: {}", queryStr);

        return mDbContext.fetch(queryStr);
    }

    private void writeTranscriptProteinData(final String outputFile)
    {
        LNX_LOGGER.info("caching protein data to {}", outputFile);

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
            LNX_LOGGER.error("error writing Ensembl trans-protein data file: {}", e.toString());
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

}
