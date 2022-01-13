package com.hartwig.hmftools.sage.panel;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.panel.RegionData.validate;
import static com.hartwig.hmftools.sage.panel.RegionType.CODING;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// Generate the panel BED file as targets for panel sequencing

public class GenerateBedRegions
{
    private final RefGenomeVersion mRefGenVersion;
    private final EnsemblDataCache mEnsemblDataCache;
    
    private final List<GeneData> mCodingGenes;
    private final List<RegionData> mSpecificRegions;
    private final List<String> mComparisonFiles;
    private final List<Integer> mTranscriptValidTSLs;

    private final Map<String,List<RegionData>> mCombinedRegions; // combined, non-overlapping regions
    private final String mSourceDir;
    private final String mOutputFile;

    private static final String SPECIFIC_REGIONS_FILE = "specific_regions_file";
    private static final String CODING_GENE_FILE = "coding_genes_file";
    private static final String TRANS_TSL_FILE = "trans_tsl_file";
    private static final String COMPARISON_BED_FILES = "comparison_bed_files";
    private static final String SOURCE_DIR = "source_dir";
    private static final String OUTPUT_FILE = "output_file";

    private static final Logger LOGGER = LogManager.getLogger(GenerateBedRegions.class);

    public GenerateBedRegions(final CommandLine cmd)
    {
        mCodingGenes = Lists.newArrayList();
        mSpecificRegions = Lists.newArrayList();
        mComparisonFiles = Lists.newArrayList();
        mTranscriptValidTSLs = Lists.newArrayList();

        mCombinedRegions = Maps.newHashMap();

        mRefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V38.toString()));
        mEnsemblDataCache = new EnsemblDataCache(cmd, mRefGenVersion);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(true);

        mSourceDir = checkAddDirSeparator(cmd.getOptionValue(SOURCE_DIR));

        try
        {
            final List<String> geneNames = Files.readAllLines(new File(mSourceDir + cmd.getOptionValue(CODING_GENE_FILE)).toPath());
            geneNames.stream()
                    .filter(x -> !x.equals("GeneName"))
                    .forEach(x -> mCodingGenes.add(mEnsemblDataCache.getGeneDataByName(x)));
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load coding genes file({}): {}", cmd.getOptionValue(CODING_GENE_FILE), e.toString());
        }

        if(cmd.hasOption(TRANS_TSL_FILE))
        {
            try
            {
                final List<String> transData = Files.readAllLines(new File(mSourceDir + cmd.getOptionValue(TRANS_TSL_FILE)).toPath());
                transData.remove(0); // header

                transData.stream().map(x -> x.split(",")).filter(x -> x.length == 2)
                        .forEach(x -> mTranscriptValidTSLs.add(Integer.parseInt(x[0])));
            }
            catch(IOException e)
            {
                LOGGER.error("failed to load Ensembl TSL file({}): {}", cmd.getOptionValue(TRANS_TSL_FILE), e.toString());
            }
        }

        loadSpecificRegions(mSourceDir + cmd.getOptionValue(SPECIFIC_REGIONS_FILE, ""));

        final List<String> geneIds = mCodingGenes.stream().map(x -> x.GeneId).collect(Collectors.toList());
        mEnsemblDataCache.loadTranscriptData(geneIds);

        mOutputFile = mSourceDir + cmd.getOptionValue(OUTPUT_FILE);

        if(cmd.hasOption(COMPARISON_BED_FILES))
        {
            mComparisonFiles.addAll(Arrays.stream(cmd.getOptionValue(COMPARISON_BED_FILES).split(";", -1)).collect(Collectors.toList()));
        }

        // ensembl_trans_tsl.csv
    }

    public void run()
    {
        LOGGER.info("generating BED file from {} coding region genes and {} specific regions",
                mCodingGenes.size(), mSpecificRegions.size());

        // first form non-overlapping coding regions from the genes
        for(GeneData geneData : mCodingGenes)
        {
            formGeneCodingRegions(geneData);
        }

        for(RegionData region : mSpecificRegions)
        {
            integrateSpecificRegion(region);
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

        writeMissingComparisonRegions();
    }

    private void formGeneCodingRegions(final GeneData geneData)
    {
        final List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(geneData.GeneId);

        for(TranscriptData transData : transcripts)
        {
            if(transData.CodingStart == null)
                continue;

            if(!transData.IsCanonical && !mTranscriptValidTSLs.isEmpty() && !mTranscriptValidTSLs.contains(transData.TransId))
                continue;

            if(transData.BioType.equals("nonsense_mediated_decay"))
                continue;

            for(ExonData exon : transData.exons())
            {
                if(exon.End < transData.CodingStart)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break;

                int regionMin = max(exon.Start, transData.CodingStart);
                int regionMax = min(exon.End, transData.CodingEnd);

                RegionData regionData = new RegionData(geneData.GeneName, new ChrBaseRegion(geneData.Chromosome, regionMin, regionMax), CODING);
                addRegion(regionData);
            }
        }
    }

    private void addRegion(final RegionData region)
    {
        List<RegionData> regions = mCombinedRegions.get(region.Region.Chromosome);

        if(regions == null)
        {
            mCombinedRegions.put(region.Region.Chromosome, Lists.newArrayList(region));
            return;
        }

        RegionData.mergeRegion(regions, region);
    }

    private void integrateSpecificRegion(final RegionData region)
    {
        List<RegionData> regions = mCombinedRegions.get(region.Region.Chromosome);

        if(regions == null)
        {
            mCombinedRegions.put(region.Region.Chromosome, Lists.newArrayList(region));
            return;
        }

        RegionData.integrateRegion(regions, region);
    }

    private void loadSpecificRegions(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

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
                    LOGGER.error("invalid specific region data: {}", line);
                    continue;
                }
            }

            LOGGER.info("loaded {} specific region entries from file: {}", mSpecificRegions.size(), filename);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load specific region file({}): {}", filename, e.toString());
        }
    }

    private void writeBedRegions()
    {
        try
        {
            final BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            int regionId = 0;

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
                    region.setId(regionId++);

                    // BED file positions require a +1 offset
                    writer.write(String.format("%s\t%d\t%d\t%s",
                            chrStr, region.Region.start() - 1, region.Region.end(), region.idName()));
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

    private void writeMissingComparisonRegions()
    {
        if(mComparisonFiles.isEmpty())
            return;

        for(String comparisonFile : mComparisonFiles)
        {
            try
            {
                List<ChrBaseRegion> comparisonRegions = Lists.newArrayList();

                final List<String> fileContents = Files.readAllLines(new File(mSourceDir + comparisonFile).toPath());

                for(final String line : fileContents)
                {
                    final String[] items = line.split("\t", -1);
                    String chromosome = items[0];
                    int posStart = Integer.parseInt(items[1]) + 1;
                    int posEnd = Integer.parseInt(items[2]);
                    comparisonRegions.add(new ChrBaseRegion(chromosome, posStart, posEnd));
                }

                LOGGER.info("loaded {} comparison regions from file: {}", comparisonRegions.size(), comparisonFile);

                String outputFile = mOutputFile.replace(".bed", "") + "_vs_"
                        + comparisonFile.replace(".bed", ".csv");

                final BufferedWriter writer = createBufferedWriter(outputFile, false);
                writer.write("RegionName,Chromosome,PosStart,PosEnd,PercCoverage");
                writer.newLine();

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
                        List<ChrBaseRegion> overlaps =
                                comparisonRegions.stream().filter(x -> x.overlaps(region.Region)).collect(Collectors.toList());

                        if(overlaps.size() == 1 && overlaps.get(0).matches(region.Region))
                        {
                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.Region.start(), region.Region.end(), 1.0));
                            writer.newLine();
                            continue;
                        }

                        if(overlaps.isEmpty())
                        {
                            // LOGGER.info("region({}) entirely missed", region);

                            //writer.write(String.format("%s\t%d\t%d\t%s",chrStr, region.Region.start() - 1, region.Region.end(), region.name()));
                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.Region.start(), region.Region.end(), 0.0));
                            writer.newLine();
                        }
                        else
                        {
                            List<ChrBaseRegion> missedRegions = Lists.newArrayList();
                            ChrBaseRegion currentSegment = null;
                            int basesCovered = 0;

                            for(int i = region.Region.start(); i <= region.Region.end(); ++i)
                            {
                                final int position = i;

                                if(overlaps.stream().anyMatch(x -> x.containsPosition(position)))
                                {
                                    ++basesCovered;
                                }
                                else
                                {
                                    if(currentSegment == null || currentSegment.end() != i - 1)
                                    {
                                        currentSegment = new ChrBaseRegion(region.Region.Chromosome, i, i);
                                        missedRegions.add(currentSegment);
                                    }
                                    else
                                    {
                                        currentSegment.setEnd(i);
                                    }
                                }
                            }

                            // LOGGER.info("region({}) partially missed", region);
                            double coverage = basesCovered / (double) region.Region.baseLength();

                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.Region.start(), region.Region.end(), coverage));
                            writer.newLine();

                            /*
                            for(BaseRegion missedRegion : missedRegions)
                            {
                                // LOGGER.info("missed region({})", missedRegion);

                                writer.write(String.format("%s\t%d\t%d\t%s",
                                        chrStr, missedRegion.start() - 1, missedRegion.end(), region.name()));
                                writer.newLine();
                            }
                            */
                        }
                    }
                }
                writer.close();
            }
            catch(IOException e)
            {
                LOGGER.error("failed to write comparison BED: {}", e.toString());
            }
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(SOURCE_DIR, true, "Path to all input and output files");
        options.addOption(SPECIFIC_REGIONS_FILE, true, "Path to the Linx cohort SVs file");
        options.addOption(CODING_GENE_FILE, true, "External LINE data sample counts");
        options.addOption(TRANS_TSL_FILE, true, "Ensembl valid TSL transcript IDs");
        options.addOption(COMPARISON_BED_FILES, true, "Comparison BED file");
        options.addOption(OUTPUT_FILE, true, "Output BED filename");
        addEnsemblDir(options);
        addRefGenomeConfig(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

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
