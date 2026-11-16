package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneRegions.buildDriverGeneRegions;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneRegions.findInvalidDriverGenes;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateExonCoverageFilename;
import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateGeneCoverageFilename;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class GeneCoverage
{
    private final Map<String,List<ExonCoverage>> mChrGeneRegionsMap;

    private static final List<Integer> DEPTH_BUCKETS = Lists.newArrayList();
    private static final List<DepthRange> DEPTH_RANGES = Lists.newArrayList();
    public static final int MAX_DEPTH_BUCKET = 10000;
    protected static final int GENE_COVERAGE_MIN_MAP_QUALITY = 10;

    public GeneCoverage(final ConfigBuilder configBuilder)
    {
        mChrGeneRegionsMap = Maps.newHashMap();

        if(configBuilder == null || !configBuilder.hasValue(ENSEMBL_DATA_DIR) || !configBuilder.hasValue(DRIVER_GENE_PANEL))
            return;

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(configBuilder);

        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);
        ensemblDataCache.createGeneNameIdMap();

        List<DriverGene> driverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);

        List<DriverGene> invalidDriverGenes = findInvalidDriverGenes(driverGenes, ensemblDataCache);
        if(!invalidDriverGenes.isEmpty())
        {
            BT_LOGGER.error("Invalid non-Ensembl gene names: {}",
                    invalidDriverGenes.stream().map( x -> x.gene()).collect(Collectors.joining(";")));
            System.exit(1);
        }

        List<String> coverageGenes = driverGenes.stream().map(x -> x.gene()).collect(Collectors.toList());

        Map<String,List<GeneRegion>> chrGeneRegions = buildDriverGeneRegions(ensemblDataCache, coverageGenes);

        for(Map.Entry<String,List<GeneRegion>> entry : chrGeneRegions.entrySet())
        {
            List<ExonCoverage> exonCoverageList = Lists.newArrayList();
            mChrGeneRegionsMap.put(entry.getKey(), exonCoverageList);

            for(GeneRegion geneRegion : entry.getValue())
            {
                exonCoverageList.add(new ExonCoverage(geneRegion));
            }
        }
    }

    public List<ExonCoverage> getGeneRegions(final ChrBaseRegion region)
    {
        List<ExonCoverage> exonCoverageList = mChrGeneRegionsMap.get(region.Chromosome);

        if(exonCoverageList == null)
            return Collections.emptyList();

        return exonCoverageList.stream().filter(x -> region.overlaps(x)).collect(Collectors.toList());
    }

    public void writeResults(final String outputDir, final String sampleId)
    {
        if(mChrGeneRegionsMap.isEmpty())
            return;

        String geneFilename = generateGeneCoverageFilename(outputDir, sampleId);
        String exonFilename = generateExonCoverageFilename(outputDir, sampleId);

        writeExonCoverage(exonFilename);
        writeGeneCoverage(geneFilename);
    }

    private void writeExonCoverage(String exonFilename)
    {
        try
        {
            List<ExonDepth> exonDepths = Lists.newArrayList();

            for(List<ExonCoverage> exonCoverageList : mChrGeneRegionsMap.values())
            {
                for(ExonCoverage exonCoverage : exonCoverageList)
                {
                    ExonDepth exonDepth = new ExonDepth(
                            exonCoverage.GeneName, exonCoverage.Chromosome, exonCoverage.start(), exonCoverage.end(),
                            exonCoverage.ExonRank, exonCoverage.medianDepth(), exonCoverage.meanDepth(),
                            toDepthCounts(exonDepthCounts(exonCoverage)));

                    exonDepths.add(exonDepth);
                }
            }

            GeneDepthFile.writeExonCoverage(exonFilename, exonDepths);
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write exon median coverage: {}", e.toString());
        }
    }

    private void writeGeneCoverage(final String geneFilename)
    {
        try
        {
            Map<String,List<ExonCoverage>> geneExonCoverageMap = Maps.newHashMap();

            for(List<ExonCoverage> exonCoverageList : mChrGeneRegionsMap.values())
            {
                for(ExonCoverage exonCoverage : exonCoverageList)
                {
                    List<ExonCoverage> geneExons = geneExonCoverageMap.get(exonCoverage.GeneName);

                    if(geneExons == null)
                    {
                        geneExons = Lists.newArrayList();
                        geneExonCoverageMap.put(exonCoverage.GeneName, geneExons);
                    }

                    geneExons.add(exonCoverage);
                }
            }

            List<GeneDepth> geneDepths = Lists.newArrayList();

            for(Map.Entry<String,List<ExonCoverage>> entry : geneExonCoverageMap.entrySet())
            {
                String geneName = entry.getKey();
                List<ExonCoverage> exonCoverageList = entry.getValue();

                Map<DepthRange,Integer> depthCounts = toDepthCounts(baseCoverageSummary(exonCoverageList));

                int minPosition = exonCoverageList.stream().mapToInt(x -> x.start()).min().orElse(0);
                int maxPosition = exonCoverageList.stream().mapToInt(x -> x.end()).max().orElse(0);

                GeneDepth geneDepth = new GeneDepth(
                        geneName, exonCoverageList.get(0).Chromosome, minPosition, maxPosition,
                        depthCounts);

                geneDepths.add(geneDepth);
            }

            GeneDepthFile.writeGeneCoverage(geneFilename, geneDepths);
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write gene coverage: {}", e.toString());
        }
    }

    /* frequency distribution of depth in units of:
        - 1 up to 30
        - 5 up to 100
        - 50 up to 500  (8 regions)
        - 100 up to 2000 (15 regions)
        - 1000 up to 10000 (8 regions)
        - 10000+ (1 region)
    */

    static
    {
        int depthBucket = 0;
        int increment = 1;
        for(; depthBucket < 30; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 10;
        for(; depthBucket < 100; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 50;
        for(; depthBucket < 500; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 100;
        for(; depthBucket < 2000; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        increment = 1000;
        for(; depthBucket <= 10000; depthBucket += increment)
        {
            DEPTH_BUCKETS.add(depthBucket);
        }

        for(int bucket = 0; bucket < DEPTH_BUCKETS.size() - 1; ++bucket)
        {
            int depth = DEPTH_BUCKETS.get(bucket);
            int depthNext = DEPTH_BUCKETS.get(bucket + 1);
            DEPTH_RANGES.add(new DepthRange(depth, depthNext - 1));
        }

        int maxDepth = DEPTH_BUCKETS.get(DEPTH_BUCKETS.size() - 1);
        DEPTH_RANGES.add(new DepthRange(maxDepth, maxDepth));
    }

    private static Map<DepthRange,Integer> toDepthCounts(final int[] baseCounts)
    {
        Map<DepthRange,Integer> depthCounts = Maps.newLinkedHashMap();

        for(int i = 0; i < baseCounts.length; ++i)
        {
            depthCounts.put(DEPTH_RANGES.get(i), baseCounts[i]);
        }

        return depthCounts;
    }

    private static int[] exonDepthCounts(final ExonCoverage exon)
    {
        int[] exonDepth = new int[DEPTH_BUCKETS.size()];

        for(int baseDepth : exon.coverage())
        {
            exonDepth[bucket(baseDepth)]++;
        }

        return exonDepth;
    }

    private static int[] baseCoverageSummary(final List<ExonCoverage> exons)
    {
        int[] geneDepth = new int[DEPTH_BUCKETS.size()];

        for(ExonCoverage exon : exons)
        {
            int[] exonDepth = exonDepthCounts(exon);

            for(int i = 0; i < geneDepth.length; ++i)
            {
                geneDepth[i] += exonDepth[i];
            }
        }

        return geneDepth;
    }

    public static int bucket(int depth)
    {
        for(int i = 0; i < DEPTH_BUCKETS.size(); ++i)
        {
            int depthBucket = DEPTH_BUCKETS.get(i);

            if(depth <= depthBucket)
                return depth < depthBucket ? i - 1 : i;
        }

        return DEPTH_BUCKETS.size() - 1;
    }
}
