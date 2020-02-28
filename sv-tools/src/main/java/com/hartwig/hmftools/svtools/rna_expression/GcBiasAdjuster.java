package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.rna_expression.GeneBamReader.DEFAULT_MIN_MAPPING_QUALITY;

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
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class GcBiasAdjuster
{
    private final RnaExpConfig mConfig;

    private final Map<String, List<Double>> mChrRegionRatios;
    private final Map<Double,List<int[]>> mRatioDepthCounts; // GC-ratio to frequency of number of reads
    private final Map<Double, Double> mGcBiasAdjustmentFactors;

    // state for each region processed
    private long mCurrentRegionStart;
    private long mCurrentRegionEnd;
    private int mCurrentReadCount;
    private int mTotalReadCount;

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
        mTotalReadCount = 0;
    }

    public boolean enabled() { return !mConfig.GcBiasFile.isEmpty(); }

    public static int positionToIndex(long position)
    {
        return (int)(position / SEGMENT_LENGTH);
    }

    public static int positionToRegion(long position)
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

    private static final int DD_DEPTH = 0;
    private static final int DD_FREQ = 1;

    public void generateDepthCounts(final Map<String, List<EnsemblGeneData>> chrGeneMap)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        BamSlicer bamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        int bamSliceCount = 0;

        for (Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneMap.entrySet())
        {
            final String chromosome = entry.getKey();

            List<Double> ratios = mChrRegionRatios.get(chromosome);

            if (ratios == null || ratios.isEmpty())
            {
                LOGGER.error("missing chr({}) ratios", chromosome);
                return;
            }

            LOGGER.debug("measuring read depth for chromosome({})", chromosome);

            long currentRegionStart = 0;
            int geneCount = 0;

            for(EnsemblGeneData geneData : entry.getValue())
            {
                long geneStartRegion = positionToRegion(geneData.GeneStart);
                currentRegionStart = max(currentRegionStart, geneStartRegion);

                ++geneCount;

                // if(geneCount > 100)
                //    break;

                while (currentRegionStart < geneData.GeneEnd)
                {
                    int regionIndex = positionToIndex(currentRegionStart);

                    double gcRatio = ratios.get(regionIndex);

                    if(gcRatio == 0)
                    {
                        currentRegionStart += SEGMENT_LENGTH;
                        continue;
                    }

                    List<int[]> depthCounts = mRatioDepthCounts.get(gcRatio);

                    if (depthCounts == null)
                    {
                        depthCounts = Lists.newArrayList();
                        mRatioDepthCounts.put(gcRatio, depthCounts);
                    }

                    mCurrentRegionEnd = currentRegionStart + SEGMENT_LENGTH - 1;
                    mCurrentReadCount = 0;

                    List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, geneData.GeneStart, geneData.GeneEnd));

                    bamSlicer.slice(samReader, regions, this::processBamRead);
                    ++bamSliceCount;

                    if((bamSliceCount % 10000) == 0)
                    {
                        LOGGER.debug("chr({}) gene({}) GC-region({}) BAM slice count({}) readCount()",
                                chromosome, geneData.GeneName, currentRegionStart, bamSliceCount, mTotalReadCount);
                    }

                    // add depth counts in ascending order
                    int index = 0;
                    boolean found = false;
                    while (index < depthCounts.size())
                    {
                        int[] depthData = depthCounts.get(index);
                        if (mCurrentReadCount == depthData[DD_DEPTH])
                        {
                            ++depthData[DD_FREQ];
                            found = true;
                            break;
                        }
                        else if (mCurrentReadCount < depthCounts.get(index)[DD_DEPTH])
                        {
                            break;
                        }

                        ++index;
                    }

                    if(!found)
                    {
                        depthCounts.add(index, new int[]{mCurrentReadCount, 1});
                    }

                    currentRegionStart += SEGMENT_LENGTH;
                }
            }
        }

        // find the median depth count for each ratio bucket
        List<Integer> medianDepthCounts = Lists.newArrayList();
        Map<Double, Integer> gcMedianCounts = Maps.newHashMap();
        Map<Double, Double> gcAvgCounts = Maps.newHashMap();

        for (Map.Entry<Double, List<int[]>> entry : mRatioDepthCounts.entrySet())
        {
            double gcRatio = entry.getKey();
            final List<int[]> depthCounts = entry.getValue();

            int depthTotal = depthCounts.stream().mapToInt(x -> x[DD_DEPTH] * x[DD_FREQ]).sum();
            int depthReads = depthCounts.stream().mapToInt(x -> x[DD_FREQ]).sum();
            double avgDepth = depthTotal / (double)depthReads;

            int medianRead = depthReads / 2;
            int total = 0;
            int medianDepthCount = 0;

            for(int[] depthData : depthCounts)
            {
                if(total + depthData[DD_FREQ] >= medianRead)
                {
                    medianDepthCount = depthData[DD_DEPTH];
                    break;
                }

                total += depthData[DD_FREQ];
            }

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
            gcAvgCounts.put(gcRatio, avgDepth);
        }

        int medianValue = medianDepthCounts.size() / 2;
        int medianDepthCount = medianDepthCounts.get(medianValue);

        // determine the adjustment factor for each ratio bucket
        for (Map.Entry<Double, Integer> entry : gcMedianCounts.entrySet())
        {
            double gcRatio = entry.getKey();
            Integer ratioMedianDepth = entry.getValue();
            double adjustmentFactor = medianDepthCount / (double)ratioMedianDepth;
            double avgDepth = gcAvgCounts.get(gcRatio);

            LOGGER.info(String.format("GC-Bias ratio(%.2f) depth(median=%d avg=%.2f) adjFactor(%.4f) vs medianDepth(%d)",
                    gcRatio, ratioMedianDepth, avgDepth, adjustmentFactor, medianDepthCount));

            mGcBiasAdjustmentFactors.put(gcRatio, adjustmentFactor);
        }

        writeGcRatioDepthCounts();
    }

    public void processBamRead(@NotNull final SAMRecord record)
    {
        if(record.getDuplicateReadFlag())
            return;

        if(record.getStart() >= mCurrentRegionStart && record.getEnd() <= mCurrentRegionEnd)
        {
            ++mCurrentReadCount;
            ++mTotalReadCount;
        }
    }

    private void writeGcRatioDepthCounts()
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("gc_depth.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GcRatio,Depth,Frequency");
            writer.newLine();

            for (Map.Entry<Double, List<int[]>> entry : mRatioDepthCounts.entrySet())
            {
                double gcRatio = entry.getKey();
                final List<int[]> depthCounts = entry.getValue();

                for (int[] depthData : depthCounts)
                {
                    writer.write(String.format("%.2f,%d,%d", gcRatio, depthData[DD_DEPTH], depthData[DD_FREQ]));
                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write GC depth frequency data file: {}", e.toString());
        }
    }

}
