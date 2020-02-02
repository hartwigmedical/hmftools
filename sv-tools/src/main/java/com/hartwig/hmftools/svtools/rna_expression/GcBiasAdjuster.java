package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenome;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class GcBiasAdjuster
{
    private final RnaExpConfig mConfig;

    private final Map<String, List<Double>> mChrRegionRatios;
    private final Map<Double,List<Integer>> mRatioDepthCounts;
    private final Map<Double, Double> mGcBiasAdjustmentFactors;

    // state for each region processed
    private long mCurrentRegionStart;
    private long mCurrentRegionEnd;
    private int mCurrentReadCount;

    private static final double RATIO_BUCKET = 0.01;
    private static final int SEGMENT_LENGTH = 1000;

    private static final Logger LOGGER = LogManager.getLogger(GcBiasAdjuster.class);

    public GcBiasAdjuster(final RnaExpConfig config)
    {
        mConfig = config;
        mChrRegionRatios = Maps.newHashMap();
        mRatioDepthCounts = Maps.newHashMap();
        mGcBiasAdjustmentFactors = Maps.newHashMap();

        mCurrentRegionStart = 0;
        mCurrentRegionEnd = 0;
        mCurrentReadCount = 0;
    }

    public boolean enabled() { return !mConfig.GcBiasFile.isEmpty(); }

    public static int positionToIndex(long position)
    {
        return (int)(position / SEGMENT_LENGTH) * SEGMENT_LENGTH;
    }

    public double getGcRatio(final String chromosome, long position)
    {
        List<Double> ratios = mChrRegionRatios.get(chromosome);

        int index = positionToIndex(position);

        if(ratios == null || index >= ratios.size())
            return -1;

        return ratios.get(index);
    }

    public double getGcBiasAdjustmentFactor(double gcRatio)
    {
        Double adjFactor = mGcBiasAdjustmentFactors.get(gcRatio);
        return adjFactor != null ? adjFactor : 1;
    }

    private static int GC_INDEX_CHR = 0;
    private static int GC_INDEX_POS = 1;
    private static int GC_INDEX_RATIO = 2;

    public void loadData()
    {
        if (!Files.exists(Paths.get(mConfig.GcBiasFile)))
        {
            generateGcCounts();
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.GcBiasFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty patient sample IDs file({})", mConfig.GcBiasFile);
            }

            String currentChr = "";
            List<Double> ratios = null;
            while (line != null)
            {
                String[] items = line.split(",");

                if (items.length != 3)
                {
                    LOGGER.warn("invalid GC-bias record: {}", line);
                    return;
                }

                final String chromosome = items[GC_INDEX_CHR];
                long position = Long.parseLong(items[GC_INDEX_POS]);
                double ratio = roundRatio(Double.parseDouble(items[GC_INDEX_RATIO]));

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;
                    ratios = Lists.newArrayList();
                    mChrRegionRatios.put(chromosome, ratios);
                }

                ratios.add(ratio);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            LOGGER.warn("failed to load GC-bias file({}): {}", mConfig.GcBiasFile, e.toString());
        }
    }

    private void generateGcCounts()
    {
        LOGGER.info("generating GC-bias ratios");

        try
        {
            final BufferedWriter writer = createBufferedWriter(mConfig.GcBiasFile, false);

            final Map<Chromosome, Long> chromsomeLengths = RefGenome.HG19.lengths();

            for (Chromosome chr : HumanChromosome.values())
            {
                final String chromosome = chr.toString();
                long length = chromsomeLengths.get(chr);

                LOGGER.debug("generating GC-bias ratios for chromosome({})", chromosome);

                List<Double> ratios = Lists.newArrayList();

                long region = 0;
                while (region < length)
                {
                    double ratio = calcGcRatio(chromosome, region);
                    ratios.add(ratio);

                    writer.write(String.format("%s,%d,%.3f", chromosome, region, ratio));
                    writer.newLine();

                    region += SEGMENT_LENGTH;
                }

                mChrRegionRatios.put(chromosome, ratios);
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write GC-bias file: {}", e.toString());
        }
    }

    private static double roundRatio(double ratio)
    {
        return round(ratio/RATIO_BUCKET) * RATIO_BUCKET;
    }

    private double calcGcRatio(final String chromosome, long regionStart)
    {
        long regionEnd = regionStart + SEGMENT_LENGTH - 1;

        try
        {
            final String bases = mConfig.RefFastaSeqFile.getSubsequenceAt(chromosome, regionStart, regionEnd).getBaseString();
            int gcCount = 0;

            for (int i = 0; i < bases.length(); ++i)
            {
                if (bases.charAt(i) == 'C' || bases.charAt(i) == 'G')
                    ++gcCount;
            }

            double ratio = gcCount / (double) SEGMENT_LENGTH;
            return roundRatio(ratio);
        }
        catch (Exception e)
        {
            return 0;
        }
    }

    public void generateDepthCounts(final RnaBamReader bamReader, final Map<String, List<EnsemblGeneData>> chrGeneMap)
    {
        int chrCount = 0;
        int bamSliceCount = 0;

        for (Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneMap.entrySet())
        {
            final String chromosome = entry.getKey();
            ++chrCount;

            List<Double> ratios = mChrRegionRatios.get(chromosome);

            if (ratios == null || ratios.isEmpty())
            {
                LOGGER.error("missing chr({}) ratios", chromosome);
                return;
            }

            LOGGER.debug("measuring read depth for chromosome({})", chromosome);

            long mCurrentRegionStart = 0;

            for(EnsemblGeneData geneData : entry.getValue())
            {
                long geneStartRegion = (geneData.GeneStart / SEGMENT_LENGTH) * SEGMENT_LENGTH;
                mCurrentRegionStart = max(mCurrentRegionStart, geneStartRegion);

                while (mCurrentRegionStart < geneData.GeneEnd)
                {
                    int regionIndex = (int)(mCurrentRegionStart/SEGMENT_LENGTH);

                    double gcRatio = ratios.get(regionIndex);

                    if(gcRatio == 0)
                    {
                        mCurrentRegionStart += SEGMENT_LENGTH;
                        continue;
                    }

                    List<Integer> depthCounts = mRatioDepthCounts.get(gcRatio);

                    if (depthCounts == null)
                    {
                        depthCounts = Lists.newArrayList();
                        mRatioDepthCounts.put(gcRatio, depthCounts);
                    }

                    // optionally cap the number of reads for each GC bucket
                    if(depthCounts.size() >= chrCount * 10)
                    {
                        mCurrentRegionStart += SEGMENT_LENGTH;
                        continue;
                    }

                    mCurrentRegionEnd = mCurrentRegionStart + SEGMENT_LENGTH - 1;
                    mCurrentReadCount = 0;

                    bamReader.readBamCounts(GenomeRegions.create(chromosome, mCurrentRegionStart, mCurrentRegionEnd), this::processBamRead);
                    ++bamSliceCount;

                    if((bamSliceCount % 1000) == 0)
                    {
                        LOGGER.debug("GC-region BAM slice count({})", bamSliceCount);
                    }

                    // add depth counts in ascending order
                    int index = 0;
                    while (index < depthCounts.size())
                    {
                        if (mCurrentReadCount < depthCounts.get(index))
                            break;

                        ++index;
                    }

                    depthCounts.add(index, mCurrentReadCount);
                }
            }
        }

        // find the median depth count for each ratio bucket
        List<Integer> medianDepthCounts = Lists.newArrayList();
        Map<Double, Integer> gcMedianCounts = Maps.newHashMap();

        for (Map.Entry<Double, List<Integer>> entry : mRatioDepthCounts.entrySet())
        {
            double gcRatio = entry.getKey();
            final List<Integer> depthCounts = entry.getValue();

            int medianIndex = depthCounts.size() / 2;
            int medianDepthCount = depthCounts.get(medianIndex);

            // add in ascending order
            int index = 0;
            while (index < medianDepthCounts.size())
            {
                if (medianDepthCount < medianDepthCounts.get(index))
                    break;

                ++index;
            }

            medianDepthCounts.add(index, medianDepthCount);
            gcMedianCounts.put(gcRatio, medianDepthCount);
        }

        int medianValue = medianDepthCounts.size() / 2;
        int medianDepthCount = medianDepthCounts.get(medianValue);

        // determine the adjustment factor for each ratio bucket
        for (Map.Entry<Double, Integer> entry : gcMedianCounts.entrySet())
        {
            double gcRatio = entry.getKey();
            Integer ratioMedianDepth = entry.getValue();
            double adjustmentFactor = medianDepthCount / (double)ratioMedianDepth;

            LOGGER.info("GC-Bias ratio({}) depth({}) adjFactor({}) vs medianDepth({})",
                    gcRatio, ratioMedianDepth, String.format("%.4f", adjustmentFactor), medianDepthCount);

            mGcBiasAdjustmentFactors.put(gcRatio, adjustmentFactor);
        }
    }

    public void processBamRead(@NotNull final SAMRecord record)
    {
        if(record.getStart() >= mCurrentRegionStart && record.getEnd() <= mCurrentRegionEnd)
        {
            ++mCurrentReadCount;
        }
    }

}
