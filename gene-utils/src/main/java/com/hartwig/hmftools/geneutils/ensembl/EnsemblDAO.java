package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_DELIM;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_GENE_DATA_FILE;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_PROTEIN_FEATURE_DATA_FILE;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_TRANS_EXON_DATA_FILE;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBaseLength;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;
import static com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache.REF_ENSEMBL_DIR;
import static com.hartwig.hmftools.geneutils.ensembl.GenerateEnsemblDataCache.HGNC_GENE_DATA_FILE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class EnsemblDAO
{
    private DSLContext mDbContext;
    private final int mCoordSystemId;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<String,GeneData> mGeneIdDataMap; // geneId to gene data
    private final Set<Integer> mTranscriptIds;

    // reference data to guide building the cache
    private final HgncGenes mHgncGenes;
    private final Map<String,GeneData> mReferenceGeneDataById;
    private final Map<String,GeneData> mReferenceGeneDataByName;
    private final Map<String,String> mGeneIdMappingOverrides;

    // reference transcripts are loaded only to evaluate difference between Ensembl versions and v37 vs v38
    private final Map<String,List<TranscriptData>> mReferenceTranscriptMap; // keyed by geneId

    private static final int COORD_SYSTEM_V37 = 2;
    private static final int COORD_SYSTEM_V38 = 4;

    // not in HGNC but retained
    private static final List<GeneData> GENE_DATA_OVERRIDES = Lists.newArrayList(
            new GeneData("ENSG00000258414", "AL121790.1","14", POS_STRAND,
                    37564047,37579125, "q21.1"));

    // GOPC processed transcript which matches a ROS1 splice site - only in v38
    private static final List<String> TRANSCRIPT_EXCLUSIONS = Lists.newArrayList("ENST00000467125");

    private static final String DB_URL = "ensembl_db";
    private static final String DB_USER = "ensembl_user";
    private static final String DB_PASS = "ensembl_pass";

    private static final String MAPPING_FILE = "gene_id_mapping_file";

    public EnsemblDAO(final CommandLine cmd)
    {
        mRefGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, String.valueOf(V37)));
        mGeneIdDataMap = Maps.newHashMap();
        mTranscriptIds = Sets.newHashSet();

        mHgncGenes = new HgncGenes(cmd.getOptionValue(HGNC_GENE_DATA_FILE));

        mReferenceTranscriptMap = Maps.newHashMap();
        mReferenceGeneDataById = Maps.newHashMap();
        mReferenceGeneDataByName = Maps.newHashMap();
        mGeneIdMappingOverrides = Maps.newHashMap();

        if(cmd.hasOption(REF_ENSEMBL_DIR))
        {
            String ensemblDataDir = cmd.getOptionValue(REF_ENSEMBL_DIR);

            Map<String, List<GeneData>> chrGeneDataMap = Maps.newHashMap();
            List<String> restrictedGeneIds = Lists.newArrayList();
            EnsemblDataLoader.loadEnsemblGeneData(ensemblDataDir, restrictedGeneIds, chrGeneDataMap, V38);

            for(Map.Entry<String, List<GeneData>> entry : chrGeneDataMap.entrySet())
            {
                for(final GeneData geneData : entry.getValue())
                {
                    mReferenceGeneDataById.put(geneData.GeneId, geneData);
                    mReferenceGeneDataByName.put(geneData.GeneName, geneData);
                }
            }

            // not currently used
            // EnsemblDataLoader.loadTranscriptData(ensemblDataDir, mReferenceTranscriptMap, restrictedGeneIds, true, true);
        }

        loadGeneIdMappings(cmd.getOptionValue(MAPPING_FILE));

        mDbContext = createEnsemblDbConnection(cmd);

        if(mDbContext == null)
        {
            mCoordSystemId = -1;
            return;
        }

        mCoordSystemId = findCoordSystemId();
        GU_LOGGER.info("refGenome({}) using coord system Id({})", mRefGenomeVersion, mCoordSystemId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - 37 (default) or 38");
        options.addOption(DB_PASS, true, "Ensembl DB password, leave out for anonymous connection");
        options.addOption(DB_URL, true, "Ensembl DB URL");
        options.addOption(DB_USER, true, "Ensembl DB username");
        options.addOption(MAPPING_FILE, true, "Optional: mapping of v37 to v38 geneIds");
    }

    public boolean isValid() { return mDbContext != null && mCoordSystemId > 0; }

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
            System.setProperty("org.jooq.no-tips", "true");

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

    private void loadGeneIdMappings(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), ENSEMBL_DELIM);
            lines.remove(0);

            int geneId37Index = fieldsIndexMap.get("GeneId37");
            int geneId38Index = fieldsIndexMap.get("GeneId38");

            for(String line : lines)
            {
                String[] values = line.split(ENSEMBL_DELIM);
                String geneId37 = values[geneId37Index];
                String geneId38 = values[geneId38Index];

                mGeneIdMappingOverrides.put(geneId37, geneId38);
            }

            GU_LOGGER.info("loaded {} gene ID mappings from file({})", mGeneIdMappingOverrides.size(), filename);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to read gene ID mappings file({}): {}", filename, e.toString());
        }
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
        GU_LOGGER.info("retrieving gene data");

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,KaryotypeBand,Synonyms");
            writer.newLine();

            String queryFile = mRefGenomeVersion == V38 && mHgncGenes.hasData() ?
                    "ensembl_sql/ensembl_hgnc_gene_data.sql" : "ensembl_sql/ensembl_gene_data.sql";

            String queryStr = readQueryString(Resources.getResource(queryFile));
            queryStr = queryStr.replaceAll("COORD_SYSTEM", String.valueOf(mCoordSystemId));

            Result<Record> results = mDbContext.fetch(queryStr);

            GU_LOGGER.debug("gene query return {} records", results.size());

            final Map<String,List<GeneData>> chrGeneMap = Maps.newHashMap();
            Map<String,GeneData> geneNameMap = Maps.newHashMap();

            for(final Record record : results)
            {
                String geneName = (String)record.get("GeneName");
                String geneId = (String)record.get("GeneId");
                String synonyms = "";

                if(mHgncGenes.hasData())
                {
                    HgncGene hgncGene = null;

                    if(mRefGenomeVersion == V38)
                    {
                        String hgncId = (String)record.get("HgncId");

                        if(hgncId.isEmpty())
                            continue;

                        hgncGene = mHgncGenes.getByHgncId(hgncId);

                        if(hgncGene == null)
                            continue;

                        geneName = hgncGene.Symbol;
                        synonyms = hgncGene.HgncId;
                    }
                    else
                    {
                        // rely on the v38 genes to find and check gene details, so v37 will limited to those genes in v38
                        GeneData refGeneData = findReferenceGeneData(geneId, geneName);

                        if(refGeneData == null)
                        {
                            // check for a manual mapping entry
                            if(!mGeneIdMappingOverrides.containsKey(geneId))
                                continue;

                            String geneId38 = mGeneIdMappingOverrides.get(geneId);

                            refGeneData = mReferenceGeneDataById.get(geneId38);

                            if(refGeneData == null)
                            {
                                GU_LOGGER.error("manual gene mapping({} -> {}) has no reference data", geneId, geneId38);
                                continue;
                            }
                        }

                        geneName = refGeneData.GeneName;
                        synonyms = (String)record.get("Synonyms");

                        if(!refGeneData.getSynonyms().isEmpty() && !synonyms.contains(refGeneData.getSynonyms()))
                            synonyms = refGeneData.getSynonyms() + ";" + synonyms;
                    }
                }
                else
                {
                    Object entrezId = record.get("EntrezId");
                    if(entrezId != null)
                        geneName = (String)entrezId;
                }

                if(mGeneIdDataMap.containsKey(geneId))
                    continue;

                // check for duplicates by name
                GeneData existingGeneData = geneNameMap.get(geneName);

                if(existingGeneData != null)
                {
                    if(mRefGenomeVersion == V38 || existingGeneData.GeneId.equals(geneId))
                        continue;

                    // some v37 genes use an older geneId and in this case favour any which are used in v38
                    if(mReferenceGeneDataById.containsKey(geneId))
                    {
                        // replace and carry on
                        geneNameMap.remove(geneName);
                        mGeneIdDataMap.remove(existingGeneData.GeneId);
                        chrGeneMap.get(existingGeneData.Chromosome).remove(existingGeneData);
                    }
                    else
                    {
                        continue;
                    }
                }

                String chromosome = (String)record.get("Chromosome");
                UInteger geneStart = (UInteger) record.get("GeneStart");
                UInteger geneEnd = (UInteger) record.get("GeneEnd");
                Byte strand = (Byte) record.get("Strand");
                String karyotypeBand = (String)record.get("KaryotypeBand");

                GeneData geneData = new GeneData(
                        geneId, geneName, chromosome, strand, geneStart.intValue(), geneEnd.intValue(), karyotypeBand);

                mGeneIdDataMap.put(geneId, geneData);
                geneNameMap.put(geneName, geneData);

                geneData.setSynonyms(synonyms);

                addGeneData(chrGeneMap, geneData);
            }

            if(mRefGenomeVersion == V38)
            {
                for(GeneData geneData : GENE_DATA_OVERRIDES)
                {
                    addGeneData(chrGeneMap, geneData);
                    mGeneIdDataMap.put(geneData.GeneId, geneData);
                }
            }

            GU_LOGGER.info("caching {} genes to {}", mGeneIdDataMap.size(), outputFile);

            for(Map.Entry<String,List<GeneData>> entry : chrGeneMap.entrySet())
            {
                for(GeneData geneData : entry.getValue())
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

    private static void addGeneData(final Map<String,List<GeneData>> chrGeneMap, final GeneData geneData)
    {
        List<GeneData> geneList = chrGeneMap.get(geneData.Chromosome);

        if(geneList == null)
        {
            chrGeneMap.put(geneData.Chromosome, Lists.newArrayList(geneData));
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

    private GeneData findReferenceGeneData(final String geneId, final String geneName)
    {
        GeneData geneData = mReferenceGeneDataById.get(geneId);

        if(geneData != null)
            return geneData;

        return mReferenceGeneDataByName.get(geneName);
    }

    private void writeTranscriptExonData(final String outputFile)
    {
        GU_LOGGER.info("retrieving transcript & exon data");

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,CanonicalTranscriptId,Strand,TransId,TransName,BioType,TransStart,TransEnd");
            writer.write(",ExonRank,ExonStart,ExonEnd,ExonPhase,ExonEndPhase,CodingStart,CodingEnd");
            writer.newLine();

            final String queryStr = readQueryString(Resources.getResource("ensembl_sql/ensembl_transcript.sql"));

            Result<Record> results = mDbContext.fetch(queryStr);

            GU_LOGGER.debug("transcript query return {} records", results.size());

            final Map<String,List<TranscriptData>> transcriptDataMap = Maps.newHashMap();
            String currentGeneId = "";
            TranscriptData currentTrans = null;
            List<TranscriptData> currentTransDataList = null;

            for(final Record record : results)
            {
                String transName = (String)record.get("Trans");

                if(TRANSCRIPT_EXCLUSIONS.contains(transName))
                    continue;

                String geneId = (String)record.get("GeneId");

                if(!currentGeneId.equals(geneId))
                {
                    if(!mGeneIdDataMap.containsKey(geneId))
                        continue;

                    currentGeneId = geneId;
                    currentTransDataList = Lists.newArrayList();
                    transcriptDataMap.put(geneId, currentTransDataList);
                }

                UInteger transId = (UInteger) record.get("TransId");
                UInteger canTransId = (UInteger) record.get("CanonicalTranscriptId");
                Integer exonRank = (Integer) record.get("ExonRank");
                UInteger exonStart = (UInteger) record.get("ExonStart");
                UInteger exonEnd = (UInteger) record.get("ExonEnd");
                Byte exonPhase = (Byte) record.get("ExonPhase");
                Byte exonPhaseEnd = (Byte) record.get("ExonEndPhase");

                if(currentTrans == null || currentTrans.TransId != transId.intValue())
                {
                    Byte strand = (Byte) record.get("Strand");
                    UInteger transStart = (UInteger) record.get("TransStart");
                    UInteger transEnd = (UInteger) record.get("TransEnd");
                    ULong codingStart = (ULong) record.get("CodingStart");
                    ULong codingEnd = (ULong) record.get("CodingEnd");
                    String bioType = (String) record.get("BioType");

                    boolean isCanonical = transId.intValue() == canTransId.intValue();

                    currentTrans = new TranscriptData(
                            transId.intValue(), transName, geneId, isCanonical, strand, transStart.intValue(), transEnd.intValue(),
                            codingStart != null ? codingStart.intValue() : null, codingEnd != null ? codingEnd.intValue() : null,
                            bioType);

                    currentTransDataList.add(currentTrans);
                }

                ExonData exonData = new ExonData(
                        transId.intValue(), exonStart.intValue(), exonEnd.intValue(), exonRank.intValue(), exonPhase, exonPhaseEnd);

                currentTrans.exons().add(exonData);

                mTranscriptIds.add(transId.intValue());
            }

            GU_LOGGER.info("caching {} transcripts & exon data to {}", mTranscriptIds.size(), outputFile);

            for(Map.Entry<String,List<TranscriptData>> entry : transcriptDataMap.entrySet())
            {
                String geneId = entry.getKey();
                List<TranscriptData> transDataList = entry.getValue();

                TranscriptData canonicalTrans = transDataList.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);
                int canonicalTransId = canonicalTrans != null ? canonicalTrans.TransId : -1;

                /*

                if(canonicalTrans == null)
                {
                    GU_LOGGER.warn("geneId({}) canonical transcript not found from {} trans", geneId, transDataList.size());
                }
                */

                for(TranscriptData transData : transDataList)
                {
                    for(ExonData exon : transData.exons())
                    {
                        writer.write(String.format("%s,%d,%d,%d,%s,%s,%d,%d",
                                geneId, canonicalTransId, transData.Strand, transData.TransId,
                                transData.TransName, transData.BioType, transData.TransStart, transData.TransEnd));

                        writer.write(String.format(",%d,%d,%d,%d,%d,%s,%s",
                                exon.Rank, exon.Start, exon.End, exon.PhaseStart, exon.PhaseEnd,
                                transData.CodingStart != null ? transData.CodingStart : "NULL",
                                transData.CodingEnd != null ? transData.CodingEnd : "NULL"));

                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl trans-exon data file: {}", e.toString());
        }
    }

    private TranscriptData findCanonicalTranscript(final String geneId, final List<TranscriptData> transDataList)
    {
        TranscriptData canonicalTrans = transDataList.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

        if(mReferenceTranscriptMap.isEmpty())
            return canonicalTrans;

        // assigment logging: GeneId,SelectedTrans,Ref38Trans,CanonicalTrans,MatchType,CodingBases

        // otherwise search for the canonical in the reference data or a match based on its exons
        List<TranscriptData> refTransDataList = mReferenceTranscriptMap.get(geneId);

        if(refTransDataList == null)
        {
            GeneData geneData = mGeneIdDataMap.get(geneId);

            GeneData altGeneData = mReferenceGeneDataByName.get(geneData.GeneName);

            if(altGeneData != null)
                refTransDataList = mReferenceTranscriptMap.get(altGeneData.GeneId);

            if(refTransDataList == null)
            {
                GU_LOGGER.debug("{},{},NONE,{},NoRefByGeneId,{},{},{}",
                        geneId, canonicalTrans.TransName, canonicalTrans.TransName,
                        codingBaseLength(canonicalTrans), -1, codingBaseLength(canonicalTrans));
                return canonicalTrans;
            }
        }

        final TranscriptData refCanonicalTrans = refTransDataList.get(0);

        // prioritisation
        // coding base + transcriptId match
        // coding bases match + is 37 ensembl canonical
        // coding bases match
        // is 37 ensembl canonical

        List<TranscriptData> matchingTrans = transDataList.stream().filter(x -> transcriptsMatch(x, refCanonicalTrans)).collect(Collectors.toList());

        if(matchingTrans.isEmpty())
        {
            GU_LOGGER.debug("{},{},{},{},NoCodingMatch,{},{},{}",
                    geneId, canonicalTrans.TransName, refCanonicalTrans.TransName, canonicalTrans.TransName,
                    codingBaseLength(canonicalTrans), codingBaseLength(refCanonicalTrans), codingBaseLength(canonicalTrans));

            return canonicalTrans;
        }

        TranscriptData matchedTrans = matchingTrans.stream().filter(x -> x.TransName.equals(refCanonicalTrans.TransName)).findFirst().orElse(null);

        if(matchedTrans != null)
        {
            if(!matchedTrans.IsCanonical)
            {
                GU_LOGGER.debug("{},{},{},{},RefNameMatch,{},{},{}",
                        geneId, matchedTrans.TransName, refCanonicalTrans.TransName, canonicalTrans.TransName,
                        codingBaseLength(matchedTrans), codingBaseLength(refCanonicalTrans), codingBaseLength(canonicalTrans));
            }

            // don't log for no change to canonical and name

            //GU_LOGGER.debug("geneId({}) choosing trans({} matching ref name as canonical over designated canonical({})",
            //        geneId, matchedTrans.TransName, canonicalTrans != null ? canonicalTrans.TransName : "none");

            return matchedTrans;
        }

        matchedTrans = matchingTrans.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

        if(matchedTrans != null)
        {
            GU_LOGGER.debug("{},{},{},{},RefCanonicalNameChange,{},{},{}",
                    geneId, matchedTrans.TransName, refCanonicalTrans.TransName, canonicalTrans.TransName,
                    codingBaseLength(matchedTrans), codingBaseLength(refCanonicalTrans), codingBaseLength(canonicalTrans));

            return matchedTrans;
        }

        // any match
        matchedTrans = matchingTrans.get(0);

        GU_LOGGER.debug("{},{},{},{},RefOtherMatch,{},{},{}",
                geneId, matchedTrans.TransName, refCanonicalTrans.TransName, canonicalTrans.TransName,
                codingBaseLength(matchedTrans), codingBaseLength(refCanonicalTrans), codingBaseLength(canonicalTrans));

        return matchedTrans;
    }

    private static boolean transcriptsMatch(final TranscriptData transData, final TranscriptData otherTransData)
    {
        if(transData.exons().size() != otherTransData.exons().size())
            return false;

        if(codingBaseLength(transData) != codingBaseLength(otherTransData))
            return false;

        for(int i = 0; i < otherTransData.exons().size(); ++i)
        {
            final ExonData exon = transData.exons().get(i);
            final ExonData otherExon = otherTransData.exons().get(i);

            if(exon.PhaseStart != otherExon.PhaseStart || exon.PhaseEnd != otherExon.PhaseEnd)
                return false;
        }

        return true;
    }

    private void writeTranscriptProteinData(final String outputFile)
    {
        GU_LOGGER.info("retrieving protein data");

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("TranscriptId,TranslationId,ProteinFeatureId,SeqStart,SeqEnd,HitDescription");
            writer.newLine();

            final String queryStr = readQueryString(Resources.getResource("ensembl_sql/ensembl_protein.sql"));
            Result<Record> results = mDbContext.fetch(queryStr);

            GU_LOGGER.info("caching protein data to {}", outputFile);

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
