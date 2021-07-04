package com.hartwig.hmftools.geneutils.mapping;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_TRANS_SPLICE_DATA_FILE;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData.SYNONYM_DELIM;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.readQueryString;
import static com.hartwig.hmftools.geneutils.mapping.MappingType.GENE_ID;
import static com.hartwig.hmftools.geneutils.mapping.MappingType.GENE_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;

public class EnsemblGeneMapper
{
    private static final String ENSEMBL_DIR_37 = "ensembl_dir_37";
    private static final String ENSEMBL_DIR_38 = "ensembl_dir_38";
    private static final String LIFT_OVER_INFO_FILE = "lift_over_file";

    private final EnsemblDataCache mGeneCache37;
    private final EnsemblDataCache mGeneCache38;
    private final BufferedWriter mWriter;

    private final Map<String,List<LiftOverRegion>> mLiftOverRegions;

    private static final int GENE_COORD_BUFFER = 100;

    public EnsemblGeneMapper(final CommandLine cmd)
    {
        GU_LOGGER.info("running Ensembl gene mapping");

        String outputDir = parseOutputDir(cmd);
        String ensemblDir37 = cmd.getOptionValue(ENSEMBL_DIR_37);
        String ensemblDir38 = cmd.getOptionValue(ENSEMBL_DIR_38);
        String liftOverFile = cmd.getOptionValue(LIFT_OVER_INFO_FILE);

        if(outputDir == null || ensemblDir37 == null || ensemblDir38 == null)
        {
            GU_LOGGER.error("missing config");
            System.exit(1);
        }

        mGeneCache37 = new EnsemblDataCache(ensemblDir37, RefGenomeVersion.V37);
        mGeneCache37.setRequireGeneSynonyms();
        mGeneCache37.load(true);

        mGeneCache38 = new EnsemblDataCache(ensemblDir38, RefGenomeVersion.V38);
        mGeneCache38.setRequireGeneSynonyms();
        mGeneCache38.load(true);
        mGeneCache38.createGeneIdDataMap();
        mGeneCache38.createGeneNameIdMap();

        mLiftOverRegions = Maps.newHashMap();
        loadLiftOverFile(liftOverFile);

        mWriter = initialiseWriter(outputDir);
    }

    public void run()
    {
        GU_LOGGER.info("running Ensembl gene mapping");

        for(Map.Entry<String,List<EnsemblGeneData>> chrEntry : mGeneCache37.getChrGeneDataMap().entrySet())
        {
            String chromosome = chrEntry.getKey();
            List<EnsemblGeneData> geneList37 = chrEntry.getValue();

            String chr38 = enforceChrPrefix(chromosome);

            List<EnsemblGeneData> geneList38 = mGeneCache38.getChrGeneDataMap().get(chr38);
            List<LiftOverRegion> liftOverRegions = mLiftOverRegions.get(chr38);

            for(EnsemblGeneData geneData37 : geneList37)
            {
                EnsemblGeneData geneData38 = mGeneCache38.getGeneDataById(geneData37.GeneId);

                if(geneData38 != null)
                {
                    writeMappingData(geneData37, geneData38, GENE_ID);
                    continue;
                }

                geneData38 = mGeneCache38.getGeneDataByName(geneData37.GeneName);

                if(geneData38 != null)
                {
                    writeMappingData(geneData37, geneData38, GENE_NAME);
                    continue;
                }

                // search manually by Entrez IDs and then gene coords
                LiftOverRegion liftOverRegion = liftOverRegions.stream().filter(x -> x.GeneId.matches(geneData37.GeneId)).findFirst().orElse(null);

                if(liftOverRegion == null)
                {
                    writeMappingData(geneData37, null, MappingType.NONE);
                    continue;
                }

                boolean found = false;

                for(EnsemblGeneData gene38 : geneList38)
                {
                    if(matchOnSynonyms(geneData37, gene38))
                    {
                        writeMappingData(geneData37, null, MappingType.SYNONYM);
                        found = true;
                        break;
                    }

                    // if(gene38.GeneStart > liftOverRegion.GeneStart + GENE_COORD_BUFFER)
                    //    continue;

                    if(liftOverRegion != null && liftOverRegion.positionMatches(gene38.GeneStart, gene38.GeneEnd, GENE_COORD_BUFFER))
                    {
                        writeMappingData(geneData37, null, MappingType.COORDS);
                        found = true;
                        break;
                    }

                    // if(gene38.GeneStart > liftOverRegion.GeneEnd)
                    //    break;
                }

                if(!found)
                    writeMappingData(geneData37, null, MappingType.NONE);
            }
        }

        closeBufferedWriter(mWriter);
        GU_LOGGER.info("Ensembl gene mapping complete");
    }

    private boolean matchOnSynonyms(final EnsemblGeneData gene37, final EnsemblGeneData gene38)
    {
        if(gene37.getSynonyms().isEmpty() || gene38.getSynonyms().isEmpty())
            return false;

        String synString2 = gene38.getSynonyms().replaceAll("HGNC:", "");

        if(gene37.getSynonyms().equals(synString2))
            return true;

        String[] synonyms1 = gene37.getSynonyms().split(SYNONYM_DELIM, -1);
        String[] synonyms2 = synString2.split(SYNONYM_DELIM, -1);

        for(int i = 0; i < synonyms1.length; ++i)
        {
            String syn1 = synonyms1[i];

            for(int j = 0; j < synonyms2.length; ++j)
            {
                String syn2 = synonyms2[j];

                if(syn1.equals(syn2))
                    return true;
            }
        }

        return false;
    }

    private void loadLiftOverFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), "\t");
            lines.remove(0);

            int geneIdIndex = fieldsIndexMap.get("GeneId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int geneStartIndex = fieldsIndexMap.get("GeneStart");
            int geneEndIndex = fieldsIndexMap.get("GeneEnd");

            for(String line : lines)
            {
                String[] items = line.split("\t");
                String chromosome = items[chromosomeIndex];

                List<LiftOverRegion> regions = mLiftOverRegions.get(chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    mLiftOverRegions.put(chromosome, regions);
                }

                regions.add(new LiftOverRegion(
                        items[geneIdIndex], Integer.parseInt(items[geneStartIndex]), Integer.parseInt(items[geneEndIndex])));
            }

            GU_LOGGER.info("loaded {} lift-over regions from file({})",
                    mLiftOverRegions.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to read gene lift-over file({}): {}", filename, e.toString());
        }
    }

    private static BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String outputFile = outputDir + "ensembl_gene_mapping.csv";

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId37,GeneName37,MappingType,GeneId38,GeneName38,OtherInfo");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error initialising Ensembl gene mapping file: {}", e.toString());
            return null;
        }
    }

    private void writeMappingData(final EnsemblGeneData geneData37, final EnsemblGeneData geneData38, MappingType type)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s", geneData37.GeneId, geneData37.GeneName, type));

            if(geneData38 != null)
            {
                String otherInfo = "";

                if(geneData37.Strand != geneData38.Strand)
                    otherInfo = "StrandMisMatch";

                mWriter.write(String.format(",%s,%s,%s",
                        geneData38.GeneId, geneData38.GeneName, otherInfo));
            }
            else
            {
                mWriter.write(",,,");
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error writing Ensembl gene mapping file: {}", e.toString());
        }
    }

    public static void main(String[] args) throws ParseException, IOException, SQLException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        EnsemblGeneMapper geneMapper = new EnsemblGeneMapper(cmd);
        geneMapper.run();
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version (V37 or V38))");
        options.addOption(ENSEMBL_DATA_DIR, true, "Path to Ensembl data cache files");
        EnsemblDAO.addCmdLineArgs(options);
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(ENSEMBL_DIR_37, true, "Ensembl data cache dir for ref-genome v37");
        options.addOption(ENSEMBL_DIR_38, true, "Ensembl data cache dir for ref-genome v38");
        options.addOption(LIFT_OVER_INFO_FILE, true, "Unmatched v37 locations lifted-over to v38");
        return options;
    }

}
