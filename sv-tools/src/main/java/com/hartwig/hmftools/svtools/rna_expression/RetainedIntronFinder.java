package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionWithin;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import ngs.Read;

public class RetainedIntronFinder
{
    private GeneReadData mGene;
    private final List<RetainedIntron> mRetainedIntrons;
    private final BufferedWriter mWriter;

    private static final int MIN_FRAG_COUNT = 3;

    public RetainedIntronFinder(final BufferedWriter writer)
    {
        mRetainedIntrons = Lists.newArrayList();
        mWriter = writer;
        mGene = null;
    }

    public void setGeneData(final GeneReadData gene)
    {
        mGene = gene;
        mRetainedIntrons.clear();
    }

    public final List<RetainedIntron> getRetainedIntrons() { return mRetainedIntrons; }

    public long[] getPositionsRange()
    {
        long[] positionsRange = new long[SE_PAIR];

        positionsRange[SE_START] = mRetainedIntrons.stream().mapToLong(x -> x.position()).min().orElse(mGene.GeneData.GeneStart);
        positionsRange[SE_END] = mRetainedIntrons.stream().mapToLong(x -> x.position()).max().orElse(mGene.GeneData.GeneEnd);

        return positionsRange;
    }

    public void evaluateFragmentReads(final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.containsSplit() || read2.containsSplit())
            return;

        // reads must span an exon boundary without being exonic in another transcript
        long spannedPosition = 0;
        boolean spannedIsStart = false;

        final List<RegionReadData> candidateRegions = Lists.newArrayList();

        for(int i = 0; i <= 1; ++i)
        {
            ReadRecord read = (i == 0) ? read1 : read2;

            if(read.getMappedRegions().values().stream().anyMatch(x -> x != EXON_INTRON))
                continue;

            for(Map.Entry<RegionReadData,RegionMatchType> entry : read.getMappedRegions().entrySet())
            {
                RegionReadData region = entry.getKey();

                if(read.getMappedRegionCoords().stream().anyMatch(x -> positionWithin(region.start(), x[SE_START], x[SE_END])))
                {
                    spannedIsStart = true;

                    // take the outer-most region(s) if there are more than one
                    if(!candidateRegions.isEmpty())
                    {
                        if(region.start() > spannedPosition)
                            continue;

                        if(spannedPosition > region.start())
                        {
                            candidateRegions.clear();
                            spannedPosition = region.start();
                        }
                    }
                    else
                    {
                        spannedPosition = region.start();
                    }

                    candidateRegions.add(region);
                }
                else if(read.getMappedRegionCoords().stream().anyMatch(x -> positionWithin(region.end(), x[SE_START], x[SE_END])))
                {
                    spannedIsStart = false;

                    // take the outer-most region(s) if there are more than one
                    if(!candidateRegions.isEmpty())
                    {
                        if(region.end() < spannedPosition)
                            continue;

                        if(spannedPosition < region.end())
                        {
                            candidateRegions.clear();
                            spannedPosition = region.end();
                        }
                    }
                    else
                    {
                        spannedPosition = region.end();
                    }

                    candidateRegions.add(region);
                }
            }
        }

        if(candidateRegions.isEmpty())
            return;

        boolean matchedExisting = false;

        for(RetainedIntron retIntron : mRetainedIntrons)
        {
            if(retIntron.position() == spannedPosition && retIntron.isStart() == spannedIsStart)
            {
                matchedExisting = true;
                retIntron.addFragmentCount();
                break;
            }
        }

        if(matchedExisting)
            return;

        RetainedIntron retIntron = new RetainedIntron(mGene.GeneData, candidateRegions, spannedIsStart);
        retIntron.addFragmentCount();
        mRetainedIntrons.add(retIntron);
    }

    public void setPositionDepthFromRead(final List<long[]> readCoords)
    {
        long readMinPos = readCoords.get(0)[SE_START];
        long readMaxPos = readCoords.get(readCoords.size() - 1)[SE_END];

        for(RetainedIntron retIntron : mRetainedIntrons)
        {
            int position = (int) retIntron.position();

            if(!positionWithin(position, readMinPos, readMaxPos))
                continue;

            if (readCoords.stream().anyMatch(x -> positionWithin(position, x[SE_START], x[SE_END])))
            {
                retIntron.addReadDepth();
            }
        }
    }

    /*
    Average unspliced coverage of gene
    Pair intron retentions with novel 5’SS, novel 3’ SS or other retained intron evidence.
     */

    public static BufferedWriter createWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("retained_intron.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,Chromosome,Strand,Position,Type,FragCount,TotalDepth,TranscriptInfo");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            RE_LOGGER.error("failed to create retained intron writer: {}", e.toString());
            return null;
        }
    }

    public void writeRetainedIntrons()
    {
        if(mWriter != null)
        {
            writeRetainedIntrons(mWriter, mRetainedIntrons);
        }
    }

    private synchronized static void writeRetainedIntrons(final BufferedWriter writer, final List<RetainedIntron> retainedIntrons)
    {
        try
        {
            for(final RetainedIntron retIntron: retainedIntrons)
            {
                if(retIntron.getFragmentCount() < MIN_FRAG_COUNT)
                    continue;

                writer.write(String.format("%s,%s,%s,%d",
                        retIntron.GeneData.GeneId,  retIntron.GeneData.GeneName,
                        retIntron.GeneData.Chromosome, retIntron.GeneData.Strand));

                writer.write(String.format(",%d,%s,%d,%d,%s",
                        retIntron.position(), retIntron.type(), retIntron.getFragmentCount(), retIntron.getDepth(),
                        retIntron.transcriptInfo()));

                writer.newLine();
            }

        }
        catch(IOException e)
        {
            RE_LOGGER.error("failed to write retained intron file: {}", e.toString());
        }
    }

}
