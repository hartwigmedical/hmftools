package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransName;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AltSpliceJunction
{
    public final GeneReadData Gene;
    public final long[] SpliceJunction;
    public final AltSpliceJunctionType Type;
    public final List<RegionReadData> StartTranscripts;
    public final List<RegionReadData> EndTranscripts;

    public final String[] RegionContexts;

    public static final String CONTEXT_SJ = "SJ";
    public static final String CONTEXT_EXONIC = "EXONIC";
    public static final String CONTEXT_INTRONIC = "INTRONIC";
    public static final String CONTEXT_MIXED = "MIXED";

    private int mFragmentCount;
    private final int[] mPositionCounts; // counts at the start and end

    /*
    Record each novel splice junction per Gene + following fields:
        distance to nearest known splice boundary at start
        distance to nearest known splice boundary at end
        Start count per (transcript combination) category
        End count per (transcript combination) category
        Annotate skipped exons, cassette exons, cryptic splice sites
        Use to analyse AHR and APC novel splice junctions in relevant samples.
        Generate PON
        Look for retained introns

     */

    public AltSpliceJunction(
            final GeneReadData geneReadData, final long[] spliceJunction, AltSpliceJunctionType type,
            final String startContext, final String endContext)
    {
        Gene = geneReadData;
        SpliceJunction = spliceJunction;

        RegionContexts = new String[] {startContext, endContext};

        StartTranscripts = Lists.newArrayList();
        EndTranscripts = Lists.newArrayList();

        Type = type;

        mFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public boolean matches(final AltSpliceJunction other)
    {
        if(!other.Gene.GeneData.GeneId.equals(Gene.GeneData.GeneId))
            return false;

        return SpliceJunction[SE_START] == other.SpliceJunction[SE_START] && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public int[] getPositionCounts() { return mPositionCounts; }
    public void setFragmentCount(int seIndex, int count) { mPositionCounts[seIndex] = count; }

    public String startTranscriptNames() { return generateTranscriptNames( StartTranscripts); }
    public String endTranscriptNames() { return generateTranscriptNames( EndTranscripts); }

    private String generateTranscriptNames(final List<RegionReadData> regions)
    {
        List<String> transNames = Lists.newArrayList();

        for(RegionReadData region: regions)
        {
            transNames.addAll(region.getRefRegions().stream().map(x -> extractTransName(x)).collect(Collectors.toList()));
        }

        return appendStrList(transNames, ';');
    }

    public long nearestStartExon()
    {
        if(!StartTranscripts.isEmpty())
            return 0;

        long nearestStartExon = -1;

        for(RegionReadData region : Gene.getExonRegions())
        {
            if(nearestStartExon != 0)
            {
                long exonDistance = abs(region.end() - SpliceJunction[SE_START]);
                nearestStartExon = nearestStartExon > 0 ? min(nearestStartExon, exonDistance) : exonDistance;
            }
        }

        return nearestStartExon;
    }

    public long nearestEndExon()
    {
        if(!EndTranscripts.isEmpty())
            return 0;

        long nearestEndExon = -1;

        for(RegionReadData region : Gene.getExonRegions())
        {
            if(nearestEndExon != 0)
            {
                long exonDistance = abs(region.start() - SpliceJunction[SE_END]);
                nearestEndExon = nearestEndExon > 0 ? min(nearestEndExon, exonDistance) : exonDistance;
            }
        }

        return nearestEndExon;
    }

    public String getBaseContext(final IndexedFastaSequenceFile refGenome, int seIndex)
    {
        long position = SpliceJunction[seIndex];
        int startOffset = (seIndex == SE_START) ? 1 : 10;
        int endOffset = startOffset == 1 ? 10: 1;

        final String baseStr = refGenome
                .getSubsequenceAt(Gene.GeneData.Chromosome, position - startOffset, position + endOffset).getBaseString();

        return baseStr;

        /*
        String reverseStr = "";

        for(int i = baseStr.length() - 1; i >= 0; --i)
        {
            reverseStr += String.valueOf(convertBase(baseStr.charAt(i)));
        }

        return reverseStr;
        */
    }



}
