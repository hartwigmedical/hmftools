package com.hartwig.hmftools.isofox.gc;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.DataUtils.sumVector;
import static com.hartwig.hmftools.common.sigs.DataUtils.sumVectors;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GC_RATIO_BUCKET;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.gc.GcRatioCounts.calcGcRatioFromReadRegions;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.exp_rates.GeneCollectionSummaryData;

public class GcTranscriptCalculator
{
    private final IsofoxConfig mConfig;

    private final EnsemblDataCache mGeneTransCache;

    private final Map<String,GcRatioCounts> mTranscriptGcRatioCache;

    private final GcRatioCounts mTranscriptFitGcCounts;
    private final double[] mGcRatioAdjustments;

    private final BufferedWriter mWriter;

    private final double[] mTotalExpectedCounts;

    public GcTranscriptCalculator(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;
        mTranscriptGcRatioCache = Maps.newHashMap();

        mTranscriptFitGcCounts = new GcRatioCounts();
        mGcRatioAdjustments = new double[mTranscriptFitGcCounts.size()];

        if(config.ExpGcRatiosFile != null)
            loadExpectedData();

        mWriter = mConfig.WriteExpectedGcRatios ? createWriter() : null;

        mTotalExpectedCounts = mConfig.WriteExpectedGcRatios ? new double[mTranscriptFitGcCounts.size()] : null;
    }

    public void close() { closeBufferedWriter(mWriter); }

    public final double[] getGcRatioAdjustments() { return mGcRatioAdjustments; }
    public final GcRatioCounts getTranscriptFitGcCounts() { return mTranscriptFitGcCounts; }

    private static final double MAX_ADJUST_FACTOR = 100;
    private static final double MIN_ADJUST_FACTOR = 1 / MAX_ADJUST_FACTOR;

    public void calcGcRatioAdjustments(final GcRatioCounts globalGcCounts)
    {
        final double[] expectedFrequencies = mTranscriptFitGcCounts.getCounts();
        final double[] actualFrequencies = globalGcCounts.getCounts();

        double actualFrequencyTotal = globalGcCounts.getCountsTotal();
        double expectedFrequencyTotal = mTranscriptFitGcCounts.getCountsTotal();

        if(expectedFrequencyTotal == 0 || actualFrequencyTotal == 0)
        {
            ISF_LOGGER.error("invalid expected({}) or actual({}) totals", expectedFrequencyTotal, actualFrequencyTotal);
            return;
        }

        for(int i = 0; i < actualFrequencies.length; ++ i)
        {
            double actualPerc = actualFrequencies[i] / actualFrequencyTotal;
            double expectedPerc = expectedFrequencies[i] / expectedFrequencyTotal;

            if(actualPerc > 0 && expectedPerc > 0)
                mGcRatioAdjustments[i] = max(min(expectedPerc/actualPerc, MAX_ADJUST_FACTOR), MIN_ADJUST_FACTOR);
            else if(actualPerc == 0 && expectedPerc == 0)
                mGcRatioAdjustments[i] = 1.0;
            else if(actualPerc == 0)
                mGcRatioAdjustments[i] = MAX_ADJUST_FACTOR;
            else
                mGcRatioAdjustments[i] = MIN_ADJUST_FACTOR;

            ISF_LOGGER.info(String.format("ratio(%.2f) actual(%.6f) expected(%.6f) adjustment(%.3f)",
                    globalGcCounts.getRatios()[i], actualPerc, expectedPerc, mGcRatioAdjustments[i]));
        }
    }

    public void generateGcCountsFromFit(final List<GeneCollectionSummaryData> geneSummaries)
    {
        // use expected GC ratio counts and 1st-pass transcript fits to derive expected GC counts
        final double[] frequencies = mTranscriptFitGcCounts.getCounts();

        for(final GeneCollectionSummaryData geneSummary : geneSummaries)
        {
            final Map<String,Double> fitAllocations = geneSummary.getFitAllocations();

            for(Map.Entry<String,Double> entry : fitAllocations.entrySet())
            {
                final String transName = entry.getKey();

                if(transName.startsWith("ENSG")) // since genes also record an allocation
                    continue;

                double fitAlloc = entry.getValue();

                final GcRatioCounts transGcCounts = mTranscriptGcRatioCache.get(transName);

                if(transGcCounts == null)
                {
                    ISF_LOGGER.error("transcript({}) missing expected GC ratio counts from cache", transName);
                    return;
                }

                final double[] transFrequencies = transGcCounts.getCounts();

                for(int i = 0; i < frequencies.length; ++i)
                {
                    frequencies[i] += transFrequencies[i] * fitAlloc;
                }
            }
        }
    }

    public void generateExpectedTranscriptCounts(final String chromosome, final List<EnsemblGeneData> geneDataList)
    {
        ISF_LOGGER.info("chromosome({}) generating expected GC ratios for {} genes", chromosome, geneDataList.size());

        for(final EnsemblGeneData geneData : geneDataList)
        {
            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            if(transDataList == null || transDataList.isEmpty())
                continue;

            transDataList.forEach(x -> calculateTranscriptGcRatios(chromosome, x));
        }

        writeExpectedGcRatios(mWriter, String.format("CHR_%s", chromosome), mTotalExpectedCounts);
    }

    private void calculateTranscriptGcRatios(final String chromosome, final TranscriptData transData)
    {
        GcRatioCounts gcRatioCounts = new GcRatioCounts();

        int readLength = mConfig.ReadLength;

        boolean endOfTrans = false;
        for (final ExonData exon : transData.exons())
        {
            for (long startPos = exon.ExonStart; startPos <= exon.ExonEnd; ++startPos)
            {
                final List<long[]> readRegions = generateReadRegions(transData, startPos, readLength);

                if(readRegions.isEmpty())
                {
                    endOfTrans = true;
                    break;
                }

                double gcRatio = calcGcRatioFromReadRegions(mConfig.RefFastaSeqFile, chromosome, readRegions);

                gcRatioCounts.addGcRatio(gcRatio);
            }

            if (endOfTrans)
                break;
        }

        sumVectors(gcRatioCounts.getCounts(), mTotalExpectedCounts);

        writeExpectedGcRatios(mWriter, transData.TransName, gcRatioCounts.getCounts());
    }

    public List<long[]> generateReadRegions(final TranscriptData transData, long startPos, int readLength)
    {
        List<long[]> readRegions = Lists.newArrayListWithExpectedSize(10);

        // set out the fragment reads either within a single exon or spanning one or more
        int exonCount = transData.exons().size();
        final ExonData lastExon = transData.exons().get(exonCount - 1);

        if(startPos + readLength - 1 > lastExon.ExonEnd)
            return readRegions;

        int remainingReadBases = readLength;
        long nextRegionStart = startPos;

        for(int i = 0; i < exonCount; ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(nextRegionStart > exon.ExonEnd)
                continue;

            if(nextRegionStart + remainingReadBases - 1 <= exon.ExonEnd)
            {
                long regionEnd = nextRegionStart + remainingReadBases - 1;
                readRegions.add(new long[] {nextRegionStart, regionEnd});
                return readRegions;
            }

            long regionEnd = exon.ExonEnd;
            long regionLength = (int)(regionEnd - nextRegionStart + 1);
            remainingReadBases -= regionLength;
            readRegions.add(new long[] {nextRegionStart, regionEnd});

            if(i == exonCount - 1)
            {
                // ran out of transcript to allocate all of this read
                readRegions.clear();
                return readRegions;
            }

            nextRegionStart = transData.exons().get(i + 1).ExonStart;
        }

        return readRegions;
    }

    private void loadExpectedData()
    {
        if(!Files.exists(Paths.get(mConfig.ExpGcRatiosFile)))
        {
            ISF_LOGGER.error("invalid expected GC ratios file");
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.ExpGcRatiosFile));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty expected GC ratios file({})", mConfig.ExpGcRatiosFile);
                return;
            }

            GcRatioCounts gcRatioCounts = new GcRatioCounts();
            int expectedColCount = 1 + gcRatioCounts.size();

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                if (items.length != expectedColCount)
                {
                    ISF_LOGGER.error("invalid exp GC ratio data length({}) vs expected({}): {}", items.length, expectedColCount, line);
                    return;
                }

                final String transName = items[0];
                gcRatioCounts = new GcRatioCounts();
                final double[] frequencies = gcRatioCounts.getCounts();

                for(int i = 1; i < items.length; ++i)
                {
                    double counts = Double.parseDouble(items[i]);
                    frequencies[i - 1] = counts;
                }

                mTranscriptGcRatioCache.put(transName, gcRatioCounts);
            }

            ISF_LOGGER.info("loaded {} transcript expected GC ratios from file({})",
                    mTranscriptGcRatioCache.size(), mConfig.ExpGcRatiosFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load expected GC ratios file({}): {}", mConfig.ExpGcRatiosFile, e.toString());return;
        }
    }

    public BufferedWriter createWriter()
    {
        try
        {
            String outputFileName = String.format("%sread_%d_%s", mConfig
                    .OutputDir, mConfig.ReadLength, "exp_gc_ratios.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("TransName");

            GcRatioCounts tmp = new GcRatioCounts();

            for(Double gcRatio : tmp.getRatios())
            {
                writer.write(String.format(",Gcr_%.2f", gcRatio));
            }

            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected GC ratio counts file: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeExpectedGcRatios(
            final BufferedWriter writer, final String transName, final double[] counts)
    {
        if(writer == null)
            return;

        try
        {
            writer.write(String.format("%s", transName));

            // convert to percentages before writing
            double frequencyTotal = sumVector(counts);

            for(Double frequency : counts)
            {
                if(frequency == 0)
                    writer.write(",0");
                else
                    writer.write(String.format(",%.6f", frequency/frequencyTotal));
            }

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected GC ratio counts file: {}", e.toString());
        }
    }

}
