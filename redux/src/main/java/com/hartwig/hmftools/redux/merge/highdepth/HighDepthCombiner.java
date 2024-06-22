package com.hartwig.hmftools.redux.merge.highdepth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.immune.ImmuneRegions.getIgRegion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class HighDepthCombiner
{
    public static final Logger MD_LOGGER = LogManager.getLogger(HighDepthCombiner.class);

    private final List<String> mInputFiles;

    private final int mMinSampleCount;
    private final int mHighDepthThreshold;
    private final String mRepeatRegionsFile;
    private final List<ChrBaseRegion> mSpecificRegions;
    private final RefGenomeVersion mRefGenVersion;
    private final List<DriverGene> mDriverGenes;
    private final KnownFusionCache mKnownFusionCache;
    private final EnsemblDataCache mEnsemblDataCache;
    private final boolean mRemoveGeneOverlaps;

    private final Map<String, List<SampleHighDepthRegion>> mChrSampleHighDepthRegions;
    private final Map<String, List<HighDepthRegion>> mFinalRegions;

    private final BufferedWriter mWriter;

    private static final String HIGH_DEPTH_FILES = "high_depth_files";
    private static final String REPEAT_REGIONS_FILE = "repeat_regions_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_SAMPLE_COUNT = "min_sample_count";
    private static final String HIGH_DEPTH_THRESHOLD = "high_depth_threshold";
    private static final String REMOVE_GENE_OVERLAPS = "remove_gene_overlaps";

    private static final int DEFAULT_MIN_SAMPLE_COUNT = 3;
    private static final int DEFAULT_HIGH_DEPTH_THRESHOLD = 250;

    public HighDepthCombiner(final ConfigBuilder configBuilder)
    {
        List<String> sampleIds = loadSampleIdsFile(configBuilder);
        String highDepthFiles = configBuilder.getValue(HIGH_DEPTH_FILES);

        mInputFiles = Lists.newArrayList();

        for(String sampleId : sampleIds)
        {
            mInputFiles.add(convertWildcardSamplePath(highDepthFiles, sampleId));
        }

        mWriter = initialiseWriter(configBuilder.getValue(OUTPUT_FILE));
        mMinSampleCount = configBuilder.getInteger(MIN_SAMPLE_COUNT);
        mHighDepthThreshold = configBuilder.getInteger(HIGH_DEPTH_THRESHOLD);
        mRepeatRegionsFile = configBuilder.hasValue(REPEAT_REGIONS_FILE) ? configBuilder.getValue(REPEAT_REGIONS_FILE) : null;
        mRemoveGeneOverlaps = configBuilder.hasFlag(REMOVE_GENE_OVERLAPS);

        mChrSampleHighDepthRegions = Maps.newHashMap();
        mFinalRegions = Maps.newHashMap();
        mRefGenVersion = RefGenomeVersion.from(configBuilder);

        mSpecificRegions = Lists.newArrayList();

        try
        {
            mSpecificRegions.addAll(loadSpecificRegions(configBuilder));
        }
        catch(ParseException e)
        {
            MD_LOGGER.error("failed to load specific regions");
        }

        if(!mRemoveGeneOverlaps)
        {
            mKnownFusionCache = null;
            mDriverGenes = null;
            mEnsemblDataCache = null;
            return;
        }

        mKnownFusionCache = new KnownFusionCache();
        mKnownFusionCache.loadFromFile(configBuilder.getValue(KNOWN_FUSIONS_FILE));

        mDriverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);

        if(configBuilder.hasValue(ENSEMBL_DATA_DIR))
        {
            mEnsemblDataCache = new EnsemblDataCache(configBuilder);
            mEnsemblDataCache.load(true);
        }
        else
        {
            mEnsemblDataCache = null;
        }
    }

    public void run()
    {
        if(mInputFiles.isEmpty())
        {
            MD_LOGGER.error("no input files specified");
            System.exit(1);
        }

        combineSampleHighDepthRegions();

        if(mRepeatRegionsFile != null)
        {
            mergeRepeatRegions();
        }

        if(mRemoveGeneOverlaps)
        {
            checkKnownGeneOverlaps();
        }

        writeCombinedResults();

        closeBufferedWriter(mWriter);

        MD_LOGGER.info("high depth region combination complete");
    }

    void mergeRepeatRegions()
    {
        MD_LOGGER.info("merging in repeat regions");

        // load and add all repeat regions with max depth set to zero.
        Map<String, List<BaseRegion>> chrRepeatRegions = loadChrBaseRegions(mRepeatRegionsFile);
        for(Map.Entry<String, List<BaseRegion>> entry : chrRepeatRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<BaseRegion> repeatRegions = entry.getValue();

            List<HighDepthRegion> regions = mFinalRegions.get(chromosome);
            if(regions == null)
            {
                regions = Lists.newArrayList();
                mFinalRegions.put(chromosome, regions);
            }

            for(BaseRegion repeatRegion : repeatRegions)
            {
                regions.add(new HighDepthRegion(new ChrBaseRegion(chromosome, repeatRegion.start(), repeatRegion.end())));
            }
        }

        // sort and merge
        for(List<HighDepthRegion> regions : mFinalRegions.values())
        {
            if(regions.size() <= 1)
            {
                continue;
            }

            Collections.sort(regions);
            HighDepthRegion prev_region = regions.get(0);
            int i = 1;
            while(i < regions.size())
            {
                HighDepthRegion region = regions.get(i);
                if(positionsOverlap(prev_region.start(), prev_region.end(), region.start(), region.end()))
                {
                    prev_region.merge(region);
                    regions.remove(i);
                    continue;
                }

                prev_region = region;
                ++i;
            }

            if(!validateRegions(regions))
            {
                System.exit(1);
            }
        }
    }

    void combineSampleHighDepthRegions()
    {
        MD_LOGGER.info("combining {} high depth region files", mInputFiles.size());

        loadSampleRegions();
        sortMergeRegions();

        for(Map.Entry<String, List<SampleHighDepthRegion>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<SampleHighDepthRegion> regions = entry.getValue();

            if(!validateRegions(regions))
            {
                System.exit(1);
            }

            List<HighDepthRegion> finalHighDepthRegions = Lists.newArrayList();
            for(SampleHighDepthRegion region : regions)
            {
                if(region.SampleMaxDepths.size() < mMinSampleCount)
                {
                    continue;
                }

                Collections.sort(region.SampleMaxDepths);
                int maxDepth = region.SampleMaxDepths.get(region.SampleMaxDepths.size() - mMinSampleCount);
                if(maxDepth >= mHighDepthThreshold)
                {
                    HighDepthRegion highDepthRegion = new HighDepthRegion(region);
                    highDepthRegion.DepthMax = maxDepth;
                    finalHighDepthRegions.add(highDepthRegion);
                }
            }

            mFinalRegions.put(chromosome, finalHighDepthRegions);
        }
    }

    private void sortMergeRegions()
    {
        for(Map.Entry<String, List<SampleHighDepthRegion>> entry : mChrSampleHighDepthRegions.entrySet())
        {
            String chromosome = entry.getKey();
            List<SampleHighDepthRegion> regions = entry.getValue();
            MD_LOGGER.info("sorting and merging regions for chromosome({})", chromosome);

            if(regions.size() <= 1)
            {
                continue;
            }

            Collections.sort(regions);
            SampleHighDepthRegion prev_region = regions.get(0);
            int i = 1;
            while(i < regions.size())
            {
                SampleHighDepthRegion region = regions.get(i);
                if(positionsOverlap(prev_region.start(), prev_region.end(), region.start(), region.end()))
                {
                    prev_region.merge(region);
                    regions.remove(i);
                    continue;
                }

                prev_region = region;
                ++i;
            }
        }
    }

    private class SampleHighDepthRegion extends ChrBaseRegion
    {
        public List<Integer> SampleMaxDepths;

        public SampleHighDepthRegion(final ChrBaseRegion region)
        {
            super(region.Chromosome, region.start(), region.end());
            SampleMaxDepths = Lists.newArrayList();
        }

        public void merge(final SampleHighDepthRegion other)
        {
            SampleMaxDepths.addAll(other.SampleMaxDepths);
            setStart(Math.min(start(), other.start()));
            setEnd(Math.max(end(), other.end()));
        }
    }

    private void checkGeneRegion(final String geneName, final Set<String> addedGenes, boolean checkUpstream)
    {
        if(addedGenes.contains(geneName))
        {
            return;
        }

        addedGenes.add(geneName);

        GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);
        if(geneData == null)
        {
            ChrBaseRegion igRegion = getIgRegion(geneName, mRefGenVersion);

            if(igRegion == null)
            {
                return;
            }

            geneData = new GeneData(geneName, geneName, igRegion.chromosome(), POS_STRAND, igRegion.start(), igRegion.end(), "");
        }

        List<HighDepthRegion> highDepthRegions = mFinalRegions.get(geneData.Chromosome);

        if(highDepthRegions == null)
        {
            return;
        }

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

            int overlappingBases =
                    max(min(geneData.GeneEnd, highDepthRegion.end()) - max(geneData.GeneStart, highDepthRegion.start()) + 1, 0);

            String overlapType;

            if(overlappingBases == 0)
            {
                overlapType = "UPSTREAM";
            }
            else if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos))
            {
                overlapType = "WITHIN_GENE";
            }
            else if(positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                overlapType = "CONTAINS_GENE";
            }
            else
            {
                overlapType = "OVERLAPS";
            }

            MD_LOGGER.debug("OVERLAP,{},{},{},{},{},{},{},{},{}",
                    geneData.GeneName, geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd,
                    highDepthRegion.start(), highDepthRegion.end(), overlappingBases, overlapType, highDepthRegion.DepthMax);

            // truncate or remove the region
            if(positionsWithin(highDepthRegion.start(), highDepthRegion.end(), geneStartPos, geneEndPos)
                    || positionsWithin(geneStartPos, geneEndPos, highDepthRegion.start(), highDepthRegion.end()))
            {
                highDepthRegions.remove(index);
                continue;
            }

            if(highDepthRegion.start() < geneStartPos)
            {
                highDepthRegion.setEnd(geneStartPos - 1);
            }
            else
            {
                highDepthRegion.setStart(geneEndPos + 1);
            }

            ++index;
        }
    }

    private void checkKnownGeneOverlaps()
    {
        if(mEnsemblDataCache == null || (mDriverGenes.isEmpty() && !mKnownFusionCache.hasValidData()))
        {
            return;
        }

        MD_LOGGER.debug("OVERLAP,GeneName,Chromosome,GeneStart,GeneEnd,RegionStart,RegionEnd,OverlapBases,OverlapType,SampleCount,DepthMin,DepthMax");

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
        MD_LOGGER.info("writing output to {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add("Chromosome");
            header.add("PosStart");
            header.add("PosEnd");
            header.add("MaxDepth");
            writer.write(header.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to initialise writer: {}", e.toString());
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
            {
                continue;
            }

            try
            {
                for(HighDepthRegion region : highDepthRegions)
                {
                    StringJoiner regionData = new StringJoiner(TSV_DELIM);
                    regionData.add(region.Chromosome);
                    regionData.add(String.valueOf(region.start()));
                    regionData.add(String.valueOf(region.end()));
                    regionData.add(String.valueOf(region.DepthMax));
                    mWriter.write(regionData.toString());
                    mWriter.newLine();
                }
            }
            catch(IOException e)
            {
                MD_LOGGER.error(" failed to write region: {}", e.toString());
            }
        }
    }

    private <E extends ChrBaseRegion> boolean validateRegions(final List<E> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            E region = regions.get(i);
            E nextRegion = regions.get(i + 1);

            if(region.end() >= nextRegion.start())
            {
                MD_LOGGER.error("region({}) overlaps with next({})", region, nextRegion);
                return false;
            }
            else if(region.start() > nextRegion.start())
            {
                MD_LOGGER.error("region({}) after with next({})", region, nextRegion);
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
                String delim = FileDelimiters.inferFileDelimiter(filename);

                lines.remove(0);

                for(String line : lines)
                {
                    String[] values = line.split(delim, -1);

                    String chromosome = values[0];

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chromosome)))
                    {
                        continue;
                    }

                    int posStart = Integer.parseInt(values[1]);
                    int posEnd = Integer.parseInt(values[2]);
                    int maxDepth = Integer.parseInt(values[4]);

                    if(!mSpecificRegions.isEmpty() && mSpecificRegions.stream().noneMatch(x ->
                            x.Chromosome.equals(chromosome) && positionsOverlap(posStart, posEnd, x.start(), x.end())))
                    {
                        continue;
                    }

                    SampleHighDepthRegion region = new SampleHighDepthRegion(new ChrBaseRegion(chromosome, posStart, posEnd));
                    region.SampleMaxDepths.add(maxDepth);

                    List<SampleHighDepthRegion> regions = mChrSampleHighDepthRegions.get(chromosome);
                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        mChrSampleHighDepthRegions.put(chromosome, regions);
                    }

                    regions.add(region);
                    ++totalRegions;
                }
            }
            catch(IOException e)
            {
                MD_LOGGER.error("failed to read high-depth regions file: {}", e.toString());
            }
        }

        MD_LOGGER.info("loaded {} high-depth regions from {} files", totalRegions, mInputFiles.size());
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(HIGH_DEPTH_FILES, true, "High depth sample file(s), use '*' in for sampleId");
        configBuilder.addConfigItem(REPEAT_REGIONS_FILE, false, "Repeat regions bed file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        configBuilder.addInteger(MIN_SAMPLE_COUNT, "Min sample count to produce region", DEFAULT_MIN_SAMPLE_COUNT);
        configBuilder.addInteger(HIGH_DEPTH_THRESHOLD, "Threshold for flagging a region as problematic", DEFAULT_HIGH_DEPTH_THRESHOLD);
        configBuilder.addFlag(REMOVE_GENE_OVERLAPS, "Remove high depth regions that overlap driver or fusion genes");
        configBuilder.addConfigItem(REF_GENOME_VERSION, REF_GENOME_VERSION_CFG_DESC);
        addGenePanelOption(configBuilder, false);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);
        addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HighDepthCombiner highDepthCombiner = new HighDepthCombiner(configBuilder);
        highDepthCombiner.run();
    }
}
