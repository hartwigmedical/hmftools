package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransId;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.extractTransName;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionWithin;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AltSpliceJunction
{
    public final GeneReadData Gene;
    public final long[] SpliceJunction;
    public final List<RegionReadData> StartRegions;
    public final List<RegionReadData> EndRegions;

    public final String[] RegionContexts;

    private AltSpliceJunctionType mType;
    private int mFragmentCount;
    private final int[] mPositionCounts; // counts at the start and end
    private final List<String> mProcessedReads;

    public static final String CONTEXT_SJ = "SPLICE_JUNC";
    public static final String CONTEXT_EXONIC = "EXONIC";
    public static final String CONTEXT_INTRONIC = "INTRONIC";
    public static final String CONTEXT_MIXED = "MIXED";

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
            final GeneReadData geneReadData, final long[] spliceJunction, AltSpliceJunctionType type, final String[] regionContexts)
    {
        Gene = geneReadData;
        SpliceJunction = spliceJunction;

        RegionContexts = regionContexts;

        StartRegions = Lists.newArrayList();
        EndRegions = Lists.newArrayList();
        mProcessedReads = Lists.newArrayList();

        mType = type;

        mFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public boolean matches(final AltSpliceJunction other)
    {
        if(!other.Gene.GeneData.GeneId.equals(Gene.GeneData.GeneId))
            return false;

        return SpliceJunction[SE_START] == other.SpliceJunction[SE_START] && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public AltSpliceJunctionType type() { return mType; }
    public void overrideType(AltSpliceJunctionType type) { mType = type; }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void setFragmentCount(int seIndex, int count) { mPositionCounts[seIndex] = count; }

    public String startTranscriptNames() { return generateTranscriptNames(StartRegions); }
    public String endTranscriptNames() { return generateTranscriptNames(EndRegions); }

    private String generateTranscriptNames(final List<RegionReadData> regions)
    {
        List<String> transNames = Lists.newArrayList();
        List<Integer> validTransIds = candidateTransIds();

        for(RegionReadData region: regions)
        {
            transNames.addAll(region.getRefRegions().stream()
                    .filter(x -> validTransIds.contains(extractTransId(x)))
                    .map(x -> extractTransName(x)).collect(Collectors.toList()));
        }

        return appendStrList(transNames, ';');
    }

    public int calcNearestExonBoundary(int seIndex)
    {
        if((seIndex == SE_START && !StartRegions.isEmpty()) || (seIndex == SE_END && !EndRegions.isEmpty()))
            return 0;

        /* find distance to nearest splice acceptor or donor as follows:
            - 5' in intron - go back to last exon end
            - 5' in exon - go forward to exon end - record as -ve value
            - 3' in intron - go forward to next exon start
            - 3' in exon - go back to exon start, record as -ve value
        */

        long position = SpliceJunction[seIndex];
        boolean forwardStrand = Gene.GeneData.Strand == 1;
        boolean isFivePrime = (seIndex == SE_START) == forwardStrand;
        boolean isExonic = RegionContexts[seIndex].equals(CONTEXT_EXONIC) || RegionContexts[seIndex].equals(CONTEXT_MIXED);
        boolean searchForwards = (isFivePrime && isExonic) || (!isFivePrime && !isExonic);

        int nearestBoundary = 0;

        for(RegionReadData region : Gene.getExonRegions())
        {
            if(isExonic && !positionWithin(position, region.start(), region.end()))
                continue;

            int distance = 0;

            if(positionWithin(position, region.start(), region.end()))
            {
                // will be negative
                distance = (int)(searchForwards ? position - region.end() : region.start() - position);
                nearestBoundary = nearestBoundary == 0 ? distance : max(distance, nearestBoundary);
            }
            else
            {
                if(isExonic)
                    continue;

                if((searchForwards && region.start() < position) || (!searchForwards && position < region.end()))
                    continue;

                // will be positive
                distance = (int)(searchForwards ? region.start() - position : position - region.end());
                nearestBoundary = nearestBoundary == 0 ? distance : min(distance, nearestBoundary);
            }
        }

        return nearestBoundary;
    }

    public String getBaseContext(final IndexedFastaSequenceFile refGenome, int seIndex)
    {
        long position = SpliceJunction[seIndex];
        int startOffset = (seIndex == SE_START) ? 1 : 10;
        int endOffset = startOffset == 1 ? 10: 1;

        final String baseStr = refGenome.getSubsequenceAt(
                Gene.GeneData.Chromosome, position - startOffset, position + endOffset).getBaseString();

        return baseStr;
    }

    public void cullNonMatchedTranscripts(final List<Integer> validTransIds)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<RegionReadData> regions = se == SE_START ? StartRegions : EndRegions;

            int index = 0;
            while(index < regions.size())
            {
                RegionReadData region = regions.get(index);

                if(validTransIds.stream().anyMatch(x -> region.hasTransId(x)))
                {
                    ++index;
                }
                else
                {
                    regions.remove(index);
                }
            }
        }
    }

    public List<Integer> candidateTransIds()
    {
        final List<Integer> uniqueTransIds = Lists.newArrayList();

        for(RegionReadData region : StartRegions)
        {
            uniqueTransIds.addAll(region.getRefRegions().stream().map(x -> extractTransId(x))
                    .filter(x -> region.getPostRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList()));
        }

        for(RegionReadData region : EndRegions)
        {
            final List<Integer> endTransIds = region.getRefRegions().stream().map(x -> extractTransId(x))
                    .filter(x -> region.getPreRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList());

            for(Integer transId : endTransIds)
            {
                if(!uniqueTransIds.contains(transId))
                    uniqueTransIds.add(transId);
            }
        }

        return uniqueTransIds;
    }


    public boolean checkProcessedRead(final String readId)
    {
        // return true and remove if found, otherwise add the new read
        for(int i = 0; i < mProcessedReads.size(); ++i)
        {
            if(mProcessedReads.get(i).equals(readId))
            {
                mProcessedReads.remove(i);
                return true;
            }
        }

        mProcessedReads.add(readId);
        return false;
    }

    public String toString()
    {
        return String.format("%s sj(%d - %d) context(%s - %s) type(%s) frags(%d)",
                Gene.GeneData.GeneId, SpliceJunction[SE_START], SpliceJunction[SE_END],
                RegionContexts[SE_START], RegionContexts[SE_END], mType, mFragmentCount);
    }

}
