package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.targetregion.RegionData.validate;
import static com.hartwig.hmftools.geneutils.targetregion.RegionType.CODING;

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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.geneutils.common.CommonUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// Generate the panel BED file as targets for panel sequencing

public class GenerateTargetRegionsBed
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
    private final boolean mIncludeUTR;

    private static final String SPECIFIC_REGIONS_FILE = "specific_regions_file";
    private static final String CODING_GENE_FILE = "coding_genes_file";
    private static final String TRANS_TSL_FILE = "trans_tsl_file";
    private static final String COMPARISON_BED_FILES = "comparison_bed_files";
    private static final String SOURCE_DIR = "source_dir";
    private static final String OUTPUT_FILE = "output_file";
    private static final String INCLUDE_UTR = "include_utr";

    private static final Logger LOGGER = LogManager.getLogger(GenerateTargetRegionsBed.class);

    public GenerateTargetRegionsBed(final ConfigBuilder configBuilder)
    {
        mCodingGenes = Lists.newArrayList();
        mSpecificRegions = Lists.newArrayList();
        mComparisonFiles = Lists.newArrayList();
        mTranscriptValidTSLs = Lists.newArrayList();

        mCombinedRegions = Maps.newHashMap();

        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(true);

        mIncludeUTR = configBuilder.hasFlag(INCLUDE_UTR);

        mSourceDir = checkAddDirSeparator(configBuilder.getValue(SOURCE_DIR));

        String codingGeneFile = configBuilder.getValue(CODING_GENE_FILE);

        try
        {
            final List<String> geneNames = Files.readAllLines(new File(mSourceDir + codingGeneFile).toPath());
            geneNames.stream()
                    .filter(x -> !x.equals("GeneName"))
                    .forEach(x -> mCodingGenes.add(mEnsemblDataCache.getGeneDataByName(x)));
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load coding genes file({}): {}", codingGeneFile, e.toString());
        }

        if(configBuilder.hasValue(TRANS_TSL_FILE))
        {
            String transTslFile = configBuilder.getValue(TRANS_TSL_FILE);
            try
            {
                final List<String> transData = Files.readAllLines(new File(mSourceDir + transTslFile).toPath());
                transData.remove(0); // header

                transData.stream().map(x -> x.split(",")).filter(x -> x.length == 2)
                        .forEach(x -> mTranscriptValidTSLs.add(Integer.parseInt(x[0])));
            }
            catch(IOException e)
            {
                LOGGER.error("failed to load Ensembl TSL file({}): {}", transTslFile, e.toString());
            }
        }

        if(configBuilder.hasValue(SPECIFIC_REGIONS_FILE))
        {
            loadSpecificRegions(mSourceDir + configBuilder.getValue(SPECIFIC_REGIONS_FILE));
        }

        final List<String> geneIds = mCodingGenes.stream().map(x -> x.GeneId).collect(Collectors.toList());
        mEnsemblDataCache.loadTranscriptData(geneIds);

        mOutputFile = mSourceDir + configBuilder.getValue(OUTPUT_FILE);

        if(configBuilder.hasValue(COMPARISON_BED_FILES))
        {
            mComparisonFiles.addAll(Arrays.stream(configBuilder.getValue(COMPARISON_BED_FILES).split(";", -1)).collect(Collectors.toList()));
        }
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

            int startPosition = mIncludeUTR || transData.nonCoding() ? transData.TransStart : transData.CodingStart;
            int endPosition = mIncludeUTR || transData.nonCoding() ? transData.TransEnd : transData.CodingEnd;

            for(ExonData exon : transData.exons())
            {
                if(exon.End < startPosition)
                    continue;

                if(exon.Start > endPosition)
                    break;

                int regionMin = max(exon.Start, startPosition);
                int regionMax = min(exon.End, endPosition);

                RegionData regionData = new RegionData(
                        geneData.GeneName, new ChrBaseRegion(geneData.Chromosome, regionMin, regionMax), exon.Rank, CODING);
                addRegion(regionData);
            }
        }
    }

    private void addRegion(final RegionData region)
    {
        List<RegionData> regions = mCombinedRegions.get(region.Chromosome);

        if(regions == null)
        {
            mCombinedRegions.put(region.Chromosome, Lists.newArrayList(region));
            return;
        }

        RegionData.mergeRegion(regions, region);
    }

    private void integrateSpecificRegion(final RegionData region)
    {
        List<RegionData> regions = mCombinedRegions.get(region.Chromosome);

        if(regions == null)
        {
            mCombinedRegions.put(region.Chromosome, Lists.newArrayList(region));
            return;
        }

        RegionData.integrateRegion(regions, region);
    }

    private void loadSpecificRegions(final String filename)
    {
        if(filename == null)
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
                    // BED file positions require a +1 offset
                    writer.write(String.format("%s\t%d\t%d\t%s",
                            chrStr, region.start() - 1, region.end(), region.idName()));
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
                                comparisonRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

                        if(overlaps.size() == 1 && overlaps.get(0).matches(region))
                        {
                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.start(), region.end(), 1.0));
                            writer.newLine();
                            continue;
                        }

                        if(overlaps.isEmpty())
                        {
                            // LOGGER.info("region({}) entirely missed", region);

                            //writer.write(String.format("%s\t%d\t%d\t%s",chrStr, region.start() - 1, region.end(), region.name()));
                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.start(), region.end(), 0.0));
                            writer.newLine();
                        }
                        else
                        {
                            List<ChrBaseRegion> missedRegions = Lists.newArrayList();
                            ChrBaseRegion currentSegment = null;
                            int basesCovered = 0;

                            for(int i = region.start(); i <= region.end(); ++i)
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
                                        currentSegment = new ChrBaseRegion(region.Chromosome, i, i);
                                        missedRegions.add(currentSegment);
                                    }
                                    else
                                    {
                                        currentSegment.setEnd(i);
                                    }
                                }
                            }

                            // LOGGER.info("region({}) partially missed", region);
                            double coverage = basesCovered / (double) region.baseLength();

                            writer.write(String.format("%s,%s,%d,%d,%.3f",
                                    region.idName(), chrStr, region.start(), region.end(), coverage));
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(SOURCE_DIR, true, "Path to all input and output files");
        configBuilder.addConfigItem(CODING_GENE_FILE, true, "External LINE data sample counts");
        configBuilder.addConfigItem(SPECIFIC_REGIONS_FILE, "Path to the Linx cohort SVs file");
        configBuilder.addConfigItem(TRANS_TSL_FILE, "Ensembl valid TSL transcript IDs");
        configBuilder.addConfigItem(COMPARISON_BED_FILES, "Comparison BED file");
        configBuilder.addFlag(INCLUDE_UTR, "Include UTR in bed regions");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output BED filename");
        addEnsemblDir(configBuilder, true);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        ConfigUtils.addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);
        CommonUtils.logVersion();

        GenerateTargetRegionsBed generateTargetRegionsBed = new GenerateTargetRegionsBed(configBuilder);
        generateTargetRegionsBed.run();

        LOGGER.info("BED region file generation complete");
    }
}
