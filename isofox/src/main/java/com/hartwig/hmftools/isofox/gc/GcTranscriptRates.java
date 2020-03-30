package com.hartwig.hmftools.isofox.gc;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.gc.GcRatioCounts.calcGcRatio;
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
import com.hartwig.hmftools.isofox.exp_rates.CategoryCountsData;

public class GcTranscriptRates
{
    private final IsofoxConfig mConfig;

    private final EnsemblDataCache mGeneTransCache;

    private final Map<String,GcRatioCounts> mTranscriptGcRatios;

    private final BufferedWriter mWriter;

    public GcTranscriptRates(final IsofoxConfig config, final EnsemblDataCache geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;
        mWriter = writer;
        mTranscriptGcRatios = Maps.newHashMap();

        if(config.ExpGcRatiosFile != null)
            loadExpectedData();
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
            int expectedColCount = 1 + gcRatioCounts.getRatios().size();

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
                final List<Double> frequencies = gcRatioCounts.getFrequencies();

                for(int i = 1; i < items.length; ++i)
                {
                    double counts = Double.parseDouble(items[i]);
                    frequencies.set(i - 1, counts);
                }

                mTranscriptGcRatios.put(transName, gcRatioCounts);
            }

            ISF_LOGGER.info("loaded {} transcript expected GC ratios from file({})",
                    mTranscriptGcRatios.size(), mConfig.ExpGcRatiosFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load expected GC ratios file({}): {}", mConfig.ExpGcRatiosFile, e.toString());
            return;
        }
    }

    public void generateExpectedRates(final String chromosome, final List<EnsemblGeneData> geneDataList)
    {
        ISF_LOGGER.info("chromosome({}) generating expected GC ratios for {} genes", chromosome, geneDataList.size());

        for(final EnsemblGeneData geneData : geneDataList)
        {
            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            if(transDataList == null || transDataList.isEmpty())
                continue;

            transDataList.forEach(x -> calculateTranscriptGcRatios(chromosome, x));
        }
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

                double gcRatio = calcGcRatioFromReadRegions(chromosome, readRegions);

                gcRatioCounts.addGcRatio(gcRatio);
            }

            if (endOfTrans)
                break;
        }

        writeExpectedGcRatios(mWriter, transData.TransName, gcRatioCounts);
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

    private double calcGcRatioFromReadRegions(final String chromosome, final List<long[]> readRegions)
    {
        double gcRatioTotal = 0;
        int basesTotal = 0;
        for(final long[] region : readRegions)
        {
            final String bases = mConfig.RefFastaSeqFile.getSubsequenceAt(chromosome, region[SE_START], region[SE_END]).getBaseString();
            basesTotal += bases.length();
            gcRatioTotal += calcGcRatio(bases) * bases.length();
        }

        return gcRatioTotal / basesTotal;
    }

    public static BufferedWriter createWriter(final IsofoxConfig config)
    {
        try
        {
            String outputFileName = String.format("%sread_%d_%s", config
                    .OutputDir, config.ReadLength, "exp_gc_ratios.csv");

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

    public synchronized static void writeExpectedGcRatios(
            final BufferedWriter writer, final String transName, final GcRatioCounts gcRatioCounts)
    {
        if(writer == null)
            return;

        try
        {
            writer.write(String.format("%s", transName));

            for(Double frequency : gcRatioCounts.getFrequencies())
            {
                writer.write(String.format(",%.1f", frequency));
            }

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected GC ratio counts file: {}", e.toString());
        }
    }

}
