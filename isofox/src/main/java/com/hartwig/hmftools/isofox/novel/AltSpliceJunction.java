package com.hartwig.hmftools.isofox.novel;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AltSpliceJunction
{
    public final long[] SpliceJunction;
    public final List<RegionReadData> StartRegions;
    public final List<RegionReadData> EndRegions;

    public final String[] RegionContexts;

    private AltSpliceJunctionType mType;
    private int mFragmentCount;
    private final int[] mPositionCounts; // counts at the start and end
    private final List<Integer> mCandidateTransIds;

    private GeneReadData mGene;

    // calculated values
    private final String[] mTranscriptNames;
    private final String[] mBaseContext;
    private final int[] mNearestExonDistance;

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

    public AltSpliceJunction(final long[] spliceJunction, AltSpliceJunctionType type, final String[] regionContexts)
    {
        SpliceJunction = spliceJunction;

        RegionContexts = regionContexts;

        StartRegions = Lists.newArrayList();
        EndRegions = Lists.newArrayList();

        mCandidateTransIds = Lists.newArrayList();
        mGene = null;

        mType = type;

        mFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
        mTranscriptNames = new String[SE_PAIR];
        mBaseContext = new String[SE_PAIR];
        mNearestExonDistance = new int[SE_PAIR];
    }

    public boolean matches(final AltSpliceJunction other)
    {
        return SpliceJunction[SE_START] == other.SpliceJunction[SE_START] && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public AltSpliceJunctionType type() { return mType; }
    public void overrideType(AltSpliceJunctionType type) { mType = type; }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void setPositionCount(int seIndex, int count) { mPositionCounts[seIndex] = count; }
    public void addPositionCount(int seIndex) { ++mPositionCounts[seIndex]; }

    public String[] getTranscriptNames() { return mTranscriptNames; }
    public String[] getBaseContext() { return mBaseContext; }
    public int[] getNearestExonDistance() { return mNearestExonDistance; }

    public void setGene(final GeneReadData gene) { mGene = gene; }
    public final GeneReadData getGene() { return mGene; }

    public void calcSummaryData(final IndexedFastaSequenceFile RefFastaSeqFile)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<RegionReadData> regions = (se == SE_START) ? StartRegions : EndRegions;
            mTranscriptNames[se] = regions.isEmpty() ? "NONE" : generateTranscriptNames(regions);
            mNearestExonDistance[se] = calcNearestExonBoundary(se);
            mBaseContext[se] = getBaseContext(RefFastaSeqFile, se);
        }
    }

    private String generateTranscriptNames(final List<RegionReadData> regions)
    {
        List<String> transNames = Lists.newArrayList();
        List<Integer> validTransIds = candidateTransIds();

        for(RegionReadData region: regions)
        {
            transNames.addAll(region.getTransExonRefs().stream()
                    .filter(x -> validTransIds.contains(x.TransId))
                    .map(x -> x.TransName).collect(Collectors.toList()));
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
        boolean forwardStrand = mGene.GeneData.Strand == 1;
        boolean isFivePrime = (seIndex == SE_START) == forwardStrand;
        boolean isExonic = RegionContexts[seIndex].equals(CONTEXT_EXONIC) || RegionContexts[seIndex].equals(CONTEXT_MIXED);
        boolean searchForwards = (isFivePrime && isExonic) || (!isFivePrime && !isExonic);

        int nearestBoundary = 0;

        for(RegionReadData region : mGene.getExonRegions())
        {
            if(isExonic && !positionWithin(position, region.start(), region.end()))
                continue;

            int distance = 0;

            if(positionWithin(position, region.start(), region.end()))
            {
                if(!region.getTransExonRefs().stream().anyMatch(x -> mCandidateTransIds.contains(x.TransId)))
                    continue;

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
                mGene.GeneData.Chromosome, position - startOffset, position + endOffset).getBaseString();

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

        int index = 0;
        while(index < mCandidateTransIds.size())
        {
            if(validTransIds.contains(mCandidateTransIds.get(index)))
            {
                ++index;
            }
            else
            {
                mCandidateTransIds.remove(index);
            }
        }
    }

    public List<Integer> candidateTransIds() { return mCandidateTransIds; }

    public void setCandidateTranscripts(final List<RegionReadData> candidateRegions)
    {
        mCandidateTransIds.clear();

        final List<Integer> validTransIds = Lists.newArrayList();

        for(RegionReadData region : candidateRegions)
        {
            if(positionWithin(SpliceJunction[SE_START], region.start(), region.end()) && positionWithin(SpliceJunction[SE_END], region.start(), region.end()))
            {
                validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId).collect(Collectors.toList()));
                continue;
            }

            if(positionWithin(SpliceJunction[SE_START], region.start(), region.end()))
            {
                // each transcript must be present in the next region to be valid
                validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId)
                        .filter(x -> region.getPostRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList()));
            }

            if(positionWithin(SpliceJunction[SE_END], region.start(), region.end()))
            {
                validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId)
                        .filter(x -> region.getPreRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList()));
            }
        }

        for(Integer transId : validTransIds)
        {
            if(!mCandidateTransIds.contains(transId))
                mCandidateTransIds.add(transId);
        }
    }

    public String toString()
    {
        return String.format("%s sj(%d - %d) context(%s - %s) type(%s) frags(%d)",
                mGene != null ? mGene.GeneData.GeneId : "unset", SpliceJunction[SE_START], SpliceJunction[SE_END],
                RegionContexts[SE_START], RegionContexts[SE_END], mType, mFragmentCount);
    }

}
