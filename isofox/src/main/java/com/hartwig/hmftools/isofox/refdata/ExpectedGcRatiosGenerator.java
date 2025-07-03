package com.hartwig.hmftools.isofox.refdata;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVectors;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.calcGcCount;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.calcGcRatio;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.calcGcRatioFromReadRegions;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.isGC;
import static com.hartwig.hmftools.isofox.refdata.RefDataWriter.writeExpectedGcRatios;

import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;

import htsjdk.samtools.SAMException;

public class ExpectedGcRatiosGenerator implements Callable<Void>
{
    private final RefDataConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final String mChromosome;
    private final List<GeneData> mGeneDataList;
    private final BufferedWriter mWriter;

    private final double[] mTotalExpectedCounts;

    public ExpectedGcRatiosGenerator(
            final RefDataConfig config, final EnsemblDataCache geneTransCache, final String chromosome, final List<GeneData> geneDataList,
            final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mGeneDataList = geneDataList;
        mChromosome = chromosome;

        mWriter = writer;

        GcRatioCounts gcRatioCounts = new GcRatioCounts();
        mTotalExpectedCounts = new double[gcRatioCounts.size()];
    }

    @Override
    public Void call()
    {
        generateExpectedCounts();
        return null;
    }

    private void generateExpectedCounts()
    {
        ISF_LOGGER.info("chromosome({}) generating expected GC ratios for {} genes", mChromosome, mGeneDataList.size());

        int genesProcessed = 0;
        int nextLogCount = 100;

        for(final GeneData geneData : mGeneDataList)
        {
            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            generateExpectedGeneCounts(geneData);

            if(transDataList != null)
            {
                transDataList.forEach(x -> calculateTranscriptGcRatios(mChromosome, x));
            }

            ++genesProcessed;

            if(genesProcessed >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chr({}) processed {} of {} genes", mChromosome, genesProcessed, mGeneDataList.size());
            }
        }

        writeExpectedGcRatios(mWriter, format("CHR_%s", mChromosome), mTotalExpectedCounts);

        ISF_LOGGER.info("chromosome({}) GC ratio generation complete", mChromosome);
    }

    private void generateExpectedGeneCounts(final GeneData geneData)
    {
        GcRatioCounts gcRatioCounts = new GcRatioCounts();
        int readLength = mConfig.ReadLength;

        final String geneBases = getRefBaseString(geneData);

        if(geneBases.length() - 1 < readLength)
        {
            gcRatioCounts.addGcRatio(calcGcRatio(geneBases));
        }
        else
        {
            // rather than measure GC content for each shifting read, just apply diffs to from the base lost and added
            List<int[]> readRegions = Lists.newArrayListWithExpectedSize(1);
            int[] geneRegion = new int[] { 0, 0 };
            readRegions.add(geneRegion);

            final String initialBases = geneBases.substring(0, readLength);
            int gcCount = calcGcCount(initialBases);
            double baseLength = mConfig.ReadLength;
            double gcRatio = gcCount / baseLength;
            gcRatioCounts.addGcRatio(gcRatio);

            for(int startPos = 1; startPos <= geneBases.length() - readLength; ++startPos)
            {
                int prevStartPos = startPos - 1;
                int endPos = startPos + readLength - 1;
                int countAdjust = (isGC(geneBases.charAt(prevStartPos)) ? -1 : 0) + (isGC(geneBases.charAt(endPos)) ? 1 : 0);

                if(gcCount > readLength || gcCount < 0)
                {
                    ISF_LOGGER.error("gene({}) gcCount error", geneData.GeneId);
                    return;
                }

                gcCount += countAdjust;
                gcRatio = gcCount / baseLength;
                gcRatioCounts.addGcRatio(gcRatio);
            }
        }

        sumVectors(gcRatioCounts.getCounts(), mTotalExpectedCounts);

        writeExpectedGcRatios(mWriter, geneData.GeneId, gcRatioCounts.getCounts());
    }

    private String getRefBaseString(final GeneData geneData)
    {
        try
        {
            return mConfig.RefGenome.getBaseString(geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd);
        }
        catch(SAMException e)
        {
            ISF_LOGGER.warn("gene({}) bases beyond ref genome", geneData);
            return "";
        }
    }

    private void calculateTranscriptGcRatios(final String chromosome, final TranscriptData transData)
    {
        GcRatioCounts gcRatioCounts = new GcRatioCounts();

        int readLength = mConfig.ReadLength;

        boolean endOfTrans = false;
        for(final ExonData exon : transData.exons())
        {
            for(int startPos = exon.Start; startPos <= exon.End; ++startPos)
            {
                final List<int[]> readRegions = generateReadRegions(transData, startPos, readLength);

                if(readRegions.isEmpty())
                {
                    endOfTrans = true;
                    break;
                }

                double gcRatio = calcGcRatioFromReadRegions(mConfig.RefGenome, chromosome, readRegions);

                gcRatioCounts.addGcRatio(gcRatio);
            }

            if(endOfTrans)
                break;
        }

        sumVectors(gcRatioCounts.getCounts(), mTotalExpectedCounts);

        writeExpectedGcRatios(mWriter, transData.TransName, gcRatioCounts.getCounts());
    }

    public List<int[]> generateReadRegions(final TranscriptData transData, int startPos, int readLength)
    {
        List<int[]> readRegions = Lists.newArrayListWithExpectedSize(10);

        // set out the fragment reads either within a single exon or spanning one or more
        int exonCount = transData.exons().size();
        final ExonData lastExon = transData.exons().get(exonCount - 1);

        if(startPos + readLength - 1 > lastExon.End)
            return readRegions;

        int remainingReadBases = readLength;
        int nextRegionStart = startPos;

        for(int i = 0; i < exonCount; ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(nextRegionStart > exon.End)
                continue;

            if(nextRegionStart + remainingReadBases - 1 <= exon.End)
            {
                int regionEnd = nextRegionStart + remainingReadBases - 1;
                readRegions.add(new int[] { nextRegionStart, regionEnd });
                return readRegions;
            }

            int regionEnd = exon.End;
            int regionLength = regionEnd - nextRegionStart + 1;
            remainingReadBases -= regionLength;
            readRegions.add(new int[] { nextRegionStart, regionEnd });

            if(i == exonCount - 1)
            {
                // ran out of transcript to allocate all of this read
                readRegions.clear();
                return readRegions;
            }

            nextRegionStart = transData.exons().get(i + 1).Start;
        }

        return readRegions;
    }
}
