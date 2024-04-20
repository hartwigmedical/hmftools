package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.enforceChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.targetregion.RegionData.validate;
import static com.hartwig.hmftools.geneutils.targetregion.RegionType.CODING;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

// Generate the panel BED file as targets for panel sequencing

public class GenerateTargetRegionsBed
{
    private final RefGenomeVersion mRefGenVersion;
    private final EnsemblDataCache mEnsemblDataCache;
    
    private final List<GeneData> mCodingGenes;
    private final List<RegionData> mSpecificRegions;
    private final List<Integer> mTranscriptValidTSLs;

    private final Map<String,List<RegionData>> mCombinedRegions; // combined, non-overlapping regions
    private final String mSourceDir;
    private final String mOutputFile;
    private final boolean mIncludeUTR;
    private final boolean mCanonicalOnly;

    private static final String SPECIFIC_REGIONS_FILE = "specific_regions_file";
    private static final String CODING_GENE_FILE = "coding_genes_file";
    private static final String TRANS_TSL_FILE = "trans_tsl_file";
    protected static final String SOURCE_DIR = "source_dir";
    protected static final String OUTPUT_FILE = "output_file";
    private static final String INCLUDE_UTR = "include_utr";
    private static final String CANONICAL_ONLY = "canonical_only";

    public GenerateTargetRegionsBed(final ConfigBuilder configBuilder)
    {
        mCodingGenes = Lists.newArrayList();
        mSpecificRegions = Lists.newArrayList();
        mTranscriptValidTSLs = Lists.newArrayList();

        mCombinedRegions = Maps.newHashMap();

        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(true);

        mIncludeUTR = configBuilder.hasFlag(INCLUDE_UTR);
        mCanonicalOnly = configBuilder.hasFlag(CANONICAL_ONLY);

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
            GU_LOGGER.error("failed to load coding genes file({}): {}", codingGeneFile, e.toString());
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
                GU_LOGGER.error("failed to load Ensembl TSL file({}): {}", transTslFile, e.toString());
            }
        }

        if(configBuilder.hasValue(SPECIFIC_REGIONS_FILE))
        {
            loadSpecificRegions(mSourceDir + configBuilder.getValue(SPECIFIC_REGIONS_FILE));
        }

        final List<String> geneIds = mCodingGenes.stream().map(x -> x.GeneId).collect(Collectors.toList());
        mEnsemblDataCache.loadTranscriptData(geneIds);

        mOutputFile = mSourceDir + configBuilder.getValue(OUTPUT_FILE);
    }

    public void run()
    {
        GU_LOGGER.info("generating BED file from {} coding region genes and {} specific regions",
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
                GU_LOGGER.error("region combining failed");
                return;
            }
        }

        GU_LOGGER.info("write {} combined regions", mCombinedRegions.values().stream().mapToInt(x -> x.size()).sum());
        writeBedRegions();

        GU_LOGGER.info("BED region file generation complete");
    }

    private void formGeneCodingRegions(final GeneData geneData)
    {
        final List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(geneData.GeneId);

        for(TranscriptData transData : transcripts)
        {
            if(transData.CodingStart == null)
                continue;

            if(!transData.IsCanonical)
            {
                if(mCanonicalOnly)
                    continue;

                if(!mTranscriptValidTSLs.isEmpty() && !mTranscriptValidTSLs.contains(transData.TransId))
                    continue;
            }

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
                    GU_LOGGER.error("invalid specific region data: {}", line);
                }
            }

            GU_LOGGER.info("loaded {} specific region entries from file: {}", mSpecificRegions.size(), filename);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to load specific region file({}): {}", filename, e.toString());
        }
    }

    private void writeBedRegions()
    {
        try
        {
            List<String> outputLines = Lists.newArrayList();

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
                    outputLines.add(format("%s\t%d\t%d\t%s",
                            chrStr, region.start() - 1, region.end(), region.idName()));
                }
            }

            NamedBedFile.write(mOutputFile, outputLines);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write to line ref-genome bases: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(SOURCE_DIR, true, "Path to all input and output files");
        configBuilder.addPrefixedPath(CODING_GENE_FILE, true, "Panel definition BED", SOURCE_DIR);
        configBuilder.addPath(SPECIFIC_REGIONS_FILE, false,"Additional regions beyond panel definition BED");
        configBuilder.addPath(TRANS_TSL_FILE, false, "Ensembl valid TSL transcript IDs");
        configBuilder.addFlag(INCLUDE_UTR, "Include UTR in bed regions");
        configBuilder.addFlag(CANONICAL_ONLY, "Form from canonical transcripts only");

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output BED filename");
        addEnsemblDir(configBuilder, true);
        addRefGenomeVersion(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateTargetRegionsBed generateTargetRegionsBed = new GenerateTargetRegionsBed(configBuilder);
        generateTargetRegionsBed.run();
    }
}
