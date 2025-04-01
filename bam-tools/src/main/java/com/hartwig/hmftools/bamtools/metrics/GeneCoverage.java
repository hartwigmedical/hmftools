package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneRegions.buildDriverGeneRegions;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateExonMediansFilename;
import static com.hartwig.hmftools.common.metrics.GeneDepthFile.generateGeneCoverageFilename;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class GeneCoverage
{
    private final Map<String,List<ExonCoverage>> mChrGeneRegionsMap;

    private static final List<Integer> DEPTH_BUCKETS = Lists.newArrayList();
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
        String exonFilename = generateExonMediansFilename(outputDir, sampleId);

        writeExonCoverage(exonFilename);
        writeGeneCoverage(geneFilename);
    }

    private void writeExonCoverage(String exonFilename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(exonFilename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("gene").add("chromosome").add("posStart").add("posEnd").add("exonRank").add("medianDepth");
            writer.write(sj.toString());
            writer.newLine();

            for(List<ExonCoverage> exonCoverageList : mChrGeneRegionsMap.values())
            {
                for(ExonCoverage exonCoverage : exonCoverageList)
                {
                    StringJoiner data = new StringJoiner(TSV_DELIM);
                    data.add(exonCoverage.GeneName);
                    data.add(exonCoverage.Chromosome);
                    data.add(String.valueOf(exonCoverage.start()));
                    data.add(String.valueOf(exonCoverage.end()));
                    data.add(String.valueOf(exonCoverage.ExonRank));
                    data.add(format("%.0f", exonCoverage.medianDepth()));
                    writer.write(data.toString());
                    writer.newLine();
                }
            }

            writer.close();
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

                int[] depthCounts = baseCoverageSummary(exonCoverageList);

                int minPosition = exonCoverageList.stream().mapToInt(x -> x.start()).min().orElse(0);
                int maxPosition = exonCoverageList.stream().mapToInt(x -> x.end()).max().orElse(0);

                GeneDepth geneDepth = new GeneDepth(
                        geneName, exonCoverageList.get(0).Chromosome, minPosition, maxPosition,
                        missedVariantLikelihood(depthCounts), depthCounts);

                geneDepths.add(geneDepth);
            }

            GeneDepthFile.write(geneFilename, geneDepths, DEPTH_BUCKETS);
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

    }

    private static int[] baseCoverageSummary(final List<ExonCoverage> exons)
    {
        int[] geneDepth = new int[DEPTH_BUCKETS.size()];

        for(ExonCoverage exon : exons)
        {
            for(int baseDepth : exon.coverage())
            {
                geneDepth[bucket(baseDepth)]++;
            }
        }

        return geneDepth;
    }

    public static int depth(int bucket)
    {
        if(bucket >= DEPTH_BUCKETS.size())
            return MAX_DEPTH_BUCKET;

        if(bucket == DEPTH_BUCKETS.size() - 1)
            return DEPTH_BUCKETS.get(DEPTH_BUCKETS.size() - 1);

        int depth = DEPTH_BUCKETS.get(bucket);
        int depthNext = DEPTH_BUCKETS.get(bucket + 1);

        if(depthNext == depth + 1)
            return depth;

        return (depth + depthNext) / 2;
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

    public static double missedVariantLikelihood(int[] baseCoverage)
    {
        int totalCoverage = Arrays.stream(baseCoverage).sum();
        double totalLikelihood = 0;

        for(int i = 0; i < baseCoverage.length; i++)
        {
            int depth = depth(i);
            int coverage = baseCoverage[i];

            if(coverage > 0)
            {
                final double proportion = 1d * coverage / totalCoverage;
                final double likelihoodOfMissing;
                if(depth == 0)
                {
                    likelihoodOfMissing = 1;
                }
                else
                {
                    final PoissonDistribution distribution = new PoissonDistribution(depth / 2d);
                    likelihoodOfMissing = distribution.cumulativeProbability(2);
                }

                totalLikelihood += proportion * likelihoodOfMissing;
            }
        }

        return totalLikelihood;
    }

}
