package com.hartwig.hmftools.svtools.bed_regions;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.svtools.bed_regions.RegionData.validate;
import static com.hartwig.hmftools.svtools.bed_regions.RegionType.CODING;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class GenerateBedRegions
{
    private final RefGenomeVersion mRefGenVersion;
    private final EnsemblDataCache mEnsemblDataCache;
    
    private final List<EnsemblGeneData> mCodingGenes;
    private final List<RegionData> mSpecificRegions;

    private final Map<String,List<RegionData>> mCombinedRegions; // combined, non-overlapping regions
    private final String mOutputDir;

    private static final String SPECIFIC_REGIONS_FILE = "specific_regions_file";
    private static final String CODING_GENE_FILE = "coding_genes_file";

    private static final Logger LOGGER = LogManager.getLogger(GenerateBedRegions.class);

    public GenerateBedRegions(final CommandLine cmd)
    {
        mCodingGenes = Lists.newArrayList();
        mSpecificRegions = Lists.newArrayList();

        mCombinedRegions = Maps.newHashMap();

        mRefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V38.toString()));
        mEnsemblDataCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), mRefGenVersion);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(true);

        try
        {
            final List<String> geneNames = Files.readAllLines(new File(cmd.getOptionValue(CODING_GENE_FILE)).toPath());
            geneNames.stream()
                    .filter(x -> !x.equals("GeneName"))
                    .forEach(x -> mCodingGenes.add(mEnsemblDataCache.getGeneDataByName(x)));
        }
        catch(IOException exception)
        {
            LOGGER.error("failed to load coding genes file({})", cmd.getOptionValue(CODING_GENE_FILE));
        }

        loadSpecificRegions(cmd.getOptionValue(SPECIFIC_REGIONS_FILE));

        final List<String> geneIds = mCodingGenes.stream().map(x -> x.GeneId).collect(Collectors.toList());
        mEnsemblDataCache.loadTranscriptData(geneIds);

        mOutputDir = parseOutputDir(cmd);
    }

    public void run()
    {
        LOGGER.info("generating BED file from {} coding region genes and {} specific regions",
                mCodingGenes.size(), mSpecificRegions.size());

        // first form non-overlapping coding regions from the genes
        for(EnsemblGeneData geneData : mCodingGenes)
        {
            formGeneCodingRegions(geneData);
        }

        for(List<RegionData> regions : mCombinedRegions.values())
        {
            if(!validate(regions))
            {
                LOGGER.error("region combining failed");
                return;
            }
        }

        LOGGER.info("write {} combined regions", mCombinedRegions.values().stream().mapToInt(x -> x.size()).sum());
        writeBedRegions();
    }

    private void formGeneCodingRegions(final EnsemblGeneData geneData)
    {
        final List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(geneData.GeneId);

        for(TranscriptData transData : transcripts)
        {
            if(transData.CodingStart == null)
                continue;

            for(ExonData exon : transData.exons())
            {
                if(exon.End < transData.CodingStart)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break;

                int regionMin = max(exon.Start, transData.CodingStart);
                int regionMax = min(exon.End, transData.CodingEnd);

                RegionData regionData = new RegionData(geneData.GeneName, new BaseRegion(geneData.Chromosome, regionMin, regionMax), CODING);
                addRegion(regionData);
            }
        }
    }

    private void addRegion(final RegionData newRegionData)
    {
        List<RegionData> regions = mCombinedRegions.get(newRegionData.Region.Chromosome);

        if(regions == null)
        {
            mCombinedRegions.put(newRegionData.Region.Chromosome, Lists.newArrayList(newRegionData));
            return;
        }

        RegionData.addRegion(regions, newRegionData);
    }

    private void loadSpecificRegions(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            final String header = fileContents.get(0);
            fileContents.remove(0);

            for(final String line : fileContents)
            {
                try
                {
                    RegionData specificRegion = RegionData.fromSpecificRegionCsv(line);

                    if(specificRegion != null)
                        mSpecificRegions.add(specificRegion);
                }
                catch(Exception e)
                {
                    LOGGER.error(" invalid specific region data: {}", line);
                    continue;
                }
            }

            LOGGER.info("loaded {} specific region entries from file: {}", mSpecificRegions.size(), filename);
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to specific region file({})", filename);
        }
    }

    private void writeBedRegions()
    {
        try
        {
            String outputFileName = mOutputDir + "bed_regions.bed";

            final BufferedWriter writer = createBufferedWriter(outputFileName, false);

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = chromosome.toString();

                if(mRefGenVersion.is38())
                    chrStr = enforceChrPrefix(chrStr);

                List<RegionData> regions = mCombinedRegions.get(chrStr);

                if(regions == null || regions.isEmpty())
                    continue;

                for(RegionData region : regions)
                {
                    writer.write(String.format("%s\t%d\t%d\t%s", chrStr, region.Region.start(), region.Region.end(), region.name()));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write to line ref-genome bases: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(SPECIFIC_REGIONS_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(CODING_GENE_FILE, true, "External LINE data sample counts");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(ENSEMBL_DATA_DIR, true, "Path to the Ensembl data cache directory");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version, 37 or 38 (default = 38)");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        GenerateBedRegions generateBedRegions = new GenerateBedRegions(cmd);
        generateBedRegions.run();

        LOGGER.info("BED region file generation complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
