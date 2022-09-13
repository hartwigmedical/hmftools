package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.tools.HighDepthConfig.HIGH_DEPTH_REGION_MAX_GAP;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.ExcludedRegions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class HighDepthCombiner
{
    private final List<String> mInputFiles;

    private final BlacklistLocations mRefefenceBlacklist;
    private final int mMinSampleCount;
    private final List<ChrBaseRegion> mSpecificRegions;
    private final RefGenomeVersion mRefGenVersion;
    private final List<DriverGene> mDriverGenes;
    private final KnownFusionCache mKnownFusionCache;
    private final EnsemblDataCache mEnsemblDataCache;
    private final boolean mRemoveGeneOverlaps;

    private final Map<String,List<List<HighDepthRegion>>> mChrSampleHighDepthRegions;
    private final Map<String,List<CombinedRegion>> mChrCombinedRegions;
    private final Map<String,List<HighDepthRegion>> mFinalRegions;

    private final BufferedWriter mWriter;

    private static final String HIGH_DEPTH_FILES = "high_depth_files";
    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BLACKLIST_FILE = "ref_blacklist_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";
    private static final String REMOVE_GENE_OVERLAPS = "remove_gene_overlaps";

    private static final int DEFAULT_MIN_SAMPLE_COUNT = 4;
    private static final int MIN_REGION_LENGTH = 11;

    public HighDepthCombiner(final CommandLine cmd)
    {
        List<String> sampleIds = loadSampleIdsFile(cmd);
        String highDepthFiles = cmd.getOptionValue(HIGH_DEPTH_FILES);

        mInputFiles = Lists.newArrayList();

        for(String sampleId : sampleIds)
        {
            mInputFiles.add(highDepthFiles.replaceAll("\\*", sampleId));
        }

        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_FILE));
        mMinSampleCount = Integer.parseInt(cmd.getOptionValue(MIN_SAMPLE_COUNT, String.valueOf(DEFAULT_MIN_SAMPLE_COUNT)));
        mRemoveGeneOverlaps = cmd.hasOption(REMOVE_GENE_OVERLAPS);

        mChrSampleHighDepthRegions = Maps.newHashMap();
        mChrCombinedRegions = Maps.newHashMap();
        mFinalRegions = Maps.newHashMap();

        mRefefenceBlacklist = new BlacklistLocations(cmd.getOptionValue(REF_BLACKLIST_FILE));
        mRefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(cmd);

        mDriverGenes = DriverGenePanelConfig.loadDriverGenes(cmd);

        if(cmd.hasOption(ENSEMBL_DATA_DIR))
        {
            mEnsemblDataCache = new EnsemblDataCache(cmd, mRefGenVersion);
            mEnsemblDataCache.load(true);
        }
        else
        {
            mEnsemblDataCache = null;
        }

        mSpecificRegions = Lists.newArrayList();

        try
        {
            mSpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            SV_LOGGER.error("failed to load specific regions");
        }
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            SV_LOGGER.error("no input files specified");
            System.exit(1);
        }

        SV_LOGGER.info("combining {} high depth region files", mInputFiles.size());

        loadSampleRegions();

        mergeRegions();

        for(Map.Entry<String,List<CombinedRegion>> entry : mChrCombinedRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<CombinedRegion> combinedRegions = entry.getValue();

            if(combinedRegions == null)
                continue;

            List<HighDepthRegion> highDepthRegions = mergeChromosomeRegions(chromosome, combinedRegions);

            if(!validateRegions(highDepthRegions))
                System.exit(1);

            mFinalRegions.put(chromosome, highDepthRegions);
        }

        checkKnownGeneOverlaps();

        writeCombinedResults();

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("high depth region combination complete");
    }

    private void mergeRegions()
    {
        for(Map.Entry<String,List<List<HighDepthRegion>>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<CombinedRegion> combinedRegions = Lists.newArrayList();
            mChrCombinedRegions.put(chromosome, combinedRegions);

            SV_LOGGER.info("merging chromosome({})", chromosome);
            int sampleIndex = 1;

            for(List<HighDepthRegion> regions : entry.getValue())
            {
                SV_LOGGER.trace("merging sample({})", sampleIndex++);

                for(HighDepthRegion region : regions)
                {
                    int index = 0;
                    boolean matched = false;

                    while(index < combinedRegions.size())
                    {
                        CombinedRegion combinedRegion = combinedRegions.get(index);
                        if(positionsOverlap(region.start(), region.end(), combinedRegion.start(), combinedRegion.end()))
                        {
                            matched = true;
                            combinedRegion.addBases(region);
                            break;
                        }
                        else if(region.end() < combinedRegion.start())
                        {
                            break;
                        }

                        ++index;
                    }

                    if(!matched)
                    {
                        CombinedRegion combinedRegion = new CombinedRegion(region);
                        combinedRegions.add(index, combinedRegion);
                    }
                    else
                    {
                        // check of this matched region now overlaps with following ones
                        CombinedRegion matchedRegion = combinedRegions.get(index);

                        int nextIndex = index + 1;
                        while(nextIndex < combinedRegions.size())
                        {
                            CombinedRegion combinedRegion = combinedRegions.get(nextIndex);

                            if(!positionsOverlap(matchedRegion.start(), matchedRegion.end(), combinedRegion.start(), combinedRegion.end()))
                                break;

                            matchedRegion.addRegion(combinedRegion);
                            combinedRegions.remove(nextIndex);
                        }
                    }
                }
            }
        }
    }

    private List<HighDepthRegion> mergeChromosomeRegions(final String chromosome, final List<CombinedRegion> combinedRegions)
    {
        List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        for(CombinedRegion region : combinedRegions)
        {
            HighDepthRegion currentRegion = null;

            for(int i = 0; i < region.Depth.size(); ++i)
            {
                PositionCount positionCount = region.Depth.get(i);

                if(positionCount.Count >= mMinSampleCount)
                {
                    if(currentRegion == null)
                    {
                        currentRegion = new HighDepthRegion(new ChrBaseRegion(chromosome, positionCount.Position, positionCount.Position));
                        currentRegion.DepthMin = positionCount.DepthMin;
                        currentRegion.DepthMax = positionCount.DepthMax;
                        currentRegion.SampleCount = positionCount.Count;
                        highDepthRegions.add(currentRegion);
                    }
                    else
                    {
                        // extend the region
                        currentRegion.setEnd(positionCount.Position);
                        currentRegion.DepthMax = max(currentRegion.DepthMax, positionCount.DepthMax);
                        currentRegion.SampleCount = max(currentRegion.SampleCount, positionCount.Count);
                    }
                }
                else
                {
                    if(currentRegion == null)
                        continue;

                    if(positionCount.Position - currentRegion.end() < HIGH_DEPTH_REGION_MAX_GAP)
                        continue;

                    // end this region
                    currentRegion = null;
                }
            }
        }

        List<BaseRegion> referenceRegions = mRefefenceBlacklist.getRegions(chromosome);

        // include the excluded region
        ChrBaseRegion excludedRegion = ExcludedRegions.getPolyGRegion(mRefGenVersion);

        if(excludedRegion.Chromosome.equals(chromosome))
        {
            referenceRegions.add(new BaseRegion(excludedRegion.start(), excludedRegion.end()));
        }

        if(referenceRegions != null)
        {
            for(BaseRegion refRegion : referenceRegions)
            {
                // merge any adjacent regions
                int index = 0;
                boolean matched = false;
                while(index < highDepthRegions.size())
                {
                    HighDepthRegion region = highDepthRegions.get(index);

                    if(region.start() > refRegion.end())
                        break;

                    if(positionsOverlap(region.start(), region.end(), refRegion.start(), refRegion.end()))
                    {
                        region.setStart(min(region.start(), refRegion.start()));
                        region.setEnd(max(region.end(), refRegion.end()));
                        matched = true;

                        // check if subsequent regions can now be merged in
                        int nextIndex = index + 1;
                        while(nextIndex < highDepthRegions.size())
                        {
                            HighDepthRegion nextRegion = highDepthRegions.get(nextIndex);

                            if(!positionsOverlap(nextRegion.start(), nextRegion.end(), refRegion.start(), refRegion.end()))
                                break;

                            highDepthRegions.remove(nextIndex);
                        }

                        break;
                    }
                    else
                    {
                        ++index;
                    }
                }

                if(!matched)
                    highDepthRegions.add(index, new HighDepthRegion(new ChrBaseRegion(chromosome, refRegion.start(), refRegion.end())));
            }
        }

        return highDepthRegions;
    }

    private class CombinedRegion
    {
        public List<PositionCount> Depth;

        public CombinedRegion(final HighDepthRegion region)
        {
            Depth = Lists.newArrayList();

            for(int pos = region.start(); pos <= region.end(); ++pos)
            {
                Depth.add(new PositionCount(pos, region.DepthMin, region.DepthMax));
            }
        }

        public int start() { return Depth.get(0).Position; }
        public int end() { return Depth.get(Depth.size() - 1).Position; }
        public int length() { return Depth.size(); }

        public void addBases(final HighDepthRegion region)
        {
            for(int pos = region.start(); pos <= region.end(); ++pos)
            {
                int existingIndex = 0;

                boolean found = false;
                while(existingIndex < Depth.size())
                {
                    PositionCount existing = Depth.get(existingIndex);

                    if(pos == existing.Position)
                    {
                        found = true;
                        ++existing.Count;
                        existing.DepthMax = max(existing.DepthMax, region.DepthMax);
                        break;
                    }

                    if(pos < existing.Position)
                        break;

                    ++existingIndex;
                }

                if(!found)
                {
                    Depth.add(existingIndex, new PositionCount(pos, region.DepthMin, region.DepthMax));
                }
            }
        }

        public void addRegion(final CombinedRegion other)
        {
            for(PositionCount otherCount : other.Depth)
            {
                int existingIndex = 0;

                boolean found = false;
                while(existingIndex < Depth.size())
                {
                    PositionCount existing = Depth.get(existingIndex);

                    if(otherCount.Position == existing.Position)
                    {
                        found = true;
                        existing.Count += otherCount.Count;
                        break;
                    }

                    if(otherCount.Position < existing.Position)
                        break;

                    ++existingIndex;
                }

                if(!found)
                {
                    Depth.add(existingIndex, otherCount);
                }
            }
        }

        public String toString() { return format("span(%d - %d) length(%d)", start(), end(), length()); }
    }

    private class PositionCount
    {
        public int Position;
        public int Count;
        public int DepthMin;
        public int DepthMax;

        public PositionCount(final int position, int depthMin, int depthMax)
        {
            Position = position;
            Count = 1;
            DepthMax = depthMax;
            DepthMin = depthMin;
        }
    }

    private void checkGeneRegion(final String geneName, final Set<String> addedGenes, boolean checkUpstream)
    {
        if(addedGenes.contains(geneName))
            return;

        addedGenes.add(geneName);

        GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);
        if(geneData == null)
        {
            if(geneName.equals("IGL"))
                geneData = new GeneData(geneName, geneName, "22", POS_STRAND, 22380474, 23265085, "");
            if(geneName.equals("IGH"))
                geneData = new GeneData(geneName, geneName, "14", POS_STRAND, 106032614, 106434161, "");
            if(geneName.equals("IGK"))
                geneData = new GeneData(geneName, geneName, "2", NEG_STRAND, 89156874, 90274235, "");
            else
                return;
        }

        List<HighDepthRegion> highDepthRegions = mFinalRegions.get(geneData.Chromosome);

        if(highDepthRegions == null)
            return;

        int geneStartPos = checkUpstream && geneData.Strand == POS_STRAND ?
                geneData.GeneStart - 10000 : geneData.GeneStart;

        int geneEndPos = checkUpstream && geneData.Strand == NEG_STRAND ? geneData.GeneEnd + 10000 : geneData.GeneEnd;

        int index = 0;
        while(index < highDepthRegions.size())
        {
            HighDepthRegion highDepthRegion = highDepthRegions.get(index);

            if(!positionsOverlap(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                ++index;
                continue;
            }

            int overlappingBases = max(min(geneData.GeneEnd, highDepthRegion.end()) - max(geneData.GeneStart, highDepthRegion.start()) + 1, 0);

            String overlapType;

            if(overlappingBases == 0)
                overlapType = "UPSTREAM";
            else if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos))
                overlapType = "WITHIN_GENE";
            else if(positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
                overlapType = "CONTAINS_GENE";
            else
                overlapType = "OVERLAPS";

            SV_LOGGER.debug("OVERLAP,{},{},{},{},{},{},{},{},{},{}",
                    geneData.GeneName, geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd,
                    highDepthRegion.start(), highDepthRegion.end(), overlappingBases, overlapType,
                    highDepthRegion.SampleCount, highDepthRegion.DepthMin, highDepthRegion.DepthMax);

            if(!mRemoveGeneOverlaps)
            {
                ++index;
                continue;
            }

            // truncate or remove the region
            if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos)
            || positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                highDepthRegions.remove(index);
                continue;
            }

            if(highDepthRegion.start() < geneStartPos)
                highDepthRegion.setEnd(geneStartPos - 1);
            else
                highDepthRegion.setStart(geneEndPos + 1);

            ++index;
        }
    }

    private void checkKnownGeneOverlaps()
    {
        if(mEnsemblDataCache == null || (mDriverGenes.isEmpty() && !mKnownFusionCache.hasValidData()))
            return;

        SV_LOGGER.debug("OVERLAP,GeneName,Chromosome,GeneStart,GeneEnd,RegionStart,RegionEnd,OverlapBases,OverlapType,SampleCount,DepthMin,DepthMax");

        // convert driver and fusion genes into regions to compare
        Set<String> addedGenes = Sets.newHashSet();

        mDriverGenes.forEach(x -> checkGeneRegion(x.gene(), addedGenes, false));

        List<KnownFusionData> knownFusionData = mKnownFusionCache.getDataByType(KNOWN_PAIR);

        if(knownFusionData != null)
        {
            knownFusionData.forEach(x -> checkGeneRegion(x.FiveGene, addedGenes, false));
            knownFusionData.forEach(x -> checkGeneRegion(x.ThreeGene, addedGenes, true));
        }

        // check the IG regions
        checkGeneRegion("IGH", addedGenes, false);
        checkGeneRegion("IGK", addedGenes, false);
        checkGeneRegion("IGL", addedGenes, false);
    }

    private BufferedWriter initialiseWriter(final String filename)
    {
        SV_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner header = new StringJoiner("\t");
            header.add("Chromosome");
            header.add("PosStart");
            header.add("PosEnd");
            header.add("SampleCount");
            header.add("DepthMin");
            header.add("DepthMax");
            writer.write(header.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to initialise writer: {}", e.toString());
        }

        return null;
    }

    private void writeCombinedResults()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());
            List<HighDepthRegion> highDepthRegions = mFinalRegions.get(chrStr);

            if(highDepthRegions == null || highDepthRegions.isEmpty())
                continue;

            try
            {
                for(HighDepthRegion region : highDepthRegions)
                {
                    if(region.baseLength() < MIN_REGION_LENGTH)
                        continue;

                    StringJoiner regionData = new StringJoiner("\t");
                    regionData.add(region.Chromosome);
                    regionData.add(String.valueOf(region.start() - 1));  // write as a BED file, so note the -1 on the start
                    regionData.add(String.valueOf(region.end()));
                    regionData.add(String.valueOf(region.SampleCount));
                    regionData.add(String.valueOf(region.DepthMin));
                    regionData.add(String.valueOf(region.DepthMax));
                    mWriter.write(regionData.toString());
                    mWriter.newLine();
                }
            }
            catch(IOException e)
            {
                SV_LOGGER.error(" failed to write region: {}", e.toString());
            }
        }
    }

    private boolean validateRegions(final List<HighDepthRegion> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            HighDepthRegion region = regions.get(i);
            HighDepthRegion nextRegion = regions.get(i + 1);

            if(region.end() >= nextRegion.start())
            {
                SV_LOGGER.error("region({}) overlaps with next({})", region, nextRegion);
                return false;
            }
            else if(region.start() > nextRegion.start())
            {
                SV_LOGGER.error("region({}) after with next({})", region, nextRegion);
                return false;
            }
        }

        return true;
    }

    private void loadSampleRegions()
    {
        int totalRegions = 0;

        for(String filename : mInputFiles)
        {
            try
            {
                List<String> lines = Files.readAllLines(Paths.get(filename));
                lines.remove(0);

                Map<String,List<HighDepthRegion>> chrRegions = Maps.newHashMap();

                for(String line : lines)
                {
                    String[] values = line.split(DELIM, -1);

                    String chromosome = values[0];

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                        continue;

                    List<HighDepthRegion> regions = chrRegions.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        chrRegions.put(chromosome, regions);
                    }

                    int posStart = Integer.parseInt(values[1]);
                    int posEnd = Integer.parseInt(values[2]);

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x ->
                            x.Chromosome.equals(chromosome) && positionsOverlap(posStart, posEnd, x.start(), x.end())))
                        continue;

                    HighDepthRegion region = new HighDepthRegion(new ChrBaseRegion(chromosome, posStart, posEnd));
                    region.DepthMin = Integer.parseInt(values[3]);
                    region.DepthMax = values.length >= 5 ? Integer.parseInt(values[4]) : region.DepthMin;
                    regions.add(region);
                    ++totalRegions;
                }

                for(Map.Entry<String,List<HighDepthRegion>> entry : chrRegions.entrySet())
                {
                    List<List<HighDepthRegion>> sampleRegions = mChrSampleHighDepthRegions.get(entry.getKey());

                    if(sampleRegions == null)
                    {
                        sampleRegions = Lists.newArrayList();
                        mChrSampleHighDepthRegions.put(entry.getKey(), sampleRegions);
                    }

                    sampleRegions.add(entry.getValue());
                }
            }
            catch(IOException e)
            {
                SV_LOGGER.error("failed to read high-depth regions file: {}", e.toString());
            }
        }

        SV_LOGGER.info("loaded {} high-depth regions from {} files", totalRegions, mInputFiles.size());
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addSampleIdFile(options);
        options.addOption(HIGH_DEPTH_FILES, true, "High depth sample file(s), use '*' in for sampleId");
        options.addOption(OUTPUT_FILE, true, "Output file");
        options.addOption(REF_BLACKLIST_FILE, true, "Reference blacklist file to include");
        options.addOption(MIN_SAMPLE_COUNT, true, "Min sample count to produce region");
        options.addOption(REMOVE_GENE_OVERLAPS, false, "Remove high depth regions that overlap driver or fusion genes");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        addGenePanelOption(false, options);
        addKnownFusionFileOption(options);
        addEnsemblDir(options);
        addOutputOptions(options);
        addLoggingOptions(options);
        addSpecificChromosomesRegionsConfig(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(cmd);
        highDepthCombiner.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
