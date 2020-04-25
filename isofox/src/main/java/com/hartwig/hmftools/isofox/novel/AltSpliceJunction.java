package com.hartwig.hmftools.isofox.novel;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.isofox.common.RnaUtils.SP_SEQ_NEG_STRAND_1;
import static com.hartwig.hmftools.isofox.common.RnaUtils.SP_SEQ_NEG_STRAND_2;
import static com.hartwig.hmftools.isofox.common.RnaUtils.SP_SEQ_POS_STRAND_1;
import static com.hartwig.hmftools.isofox.common.RnaUtils.SP_SEQ_POS_STRAND_2;
import static com.hartwig.hmftools.isofox.common.RnaUtils.formDonorAcceptorBases;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RnaUtils.setJunctionBaseContext;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionContext.EXONIC;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionContext.MIXED;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class AltSpliceJunction
{
    public final String Chromosome;
    public final int[] SpliceJunction;

    private final List<RegionReadData> mSjStartRegions; // regions which match this alt-SJ at the start
    private final List<RegionReadData> mSjEndRegions;

    public final AltSpliceJunctionContext[] RegionContexts;

    private AltSpliceJunctionType mType;
    private int mFragmentCount;
    private final int[] mPositionCounts; // counts at the start and end
    private final List<Integer> mCandidateTransIds;

    private String mGeneId; // associated gene if known or prioritised from amongst a set of candidates

    // calculated values
    private final String[] mTranscriptNames;
    private final String[] mBaseContext;
    private String mDonorAcceptorBases;
    private final int[] mNearestExonDistance;

    public AltSpliceJunction(
            final String chromosome, final int[] spliceJunction, AltSpliceJunctionType type,
            final AltSpliceJunctionContext[] regionContexts, final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions)
    {
        Chromosome = chromosome;
        SpliceJunction = spliceJunction;
        RegionContexts = regionContexts;

        mSjStartRegions = sjStartRegions;
        mSjEndRegions = sjEndRegions;

        mCandidateTransIds = Lists.newArrayList();
        mGeneId = null;

        mType = type;

        mFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
        mTranscriptNames = new String[SE_PAIR];
        mBaseContext = new String[SE_PAIR];
        mNearestExonDistance = new int[SE_PAIR];
        mDonorAcceptorBases = "";
    }

    public boolean matches(final AltSpliceJunction other)
    {
        return Chromosome.equals(other.Chromosome)
                &&  SpliceJunction[SE_START] == other.SpliceJunction[SE_START]
                && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public AltSpliceJunctionType type() { return mType; }
    public void overrideType(AltSpliceJunctionType type) { mType = type; }

    public final List<RegionReadData> getSjStartRegions() { return mSjStartRegions; }
    public final List<RegionReadData> getSjEndRegions() { return mSjEndRegions; }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public void addFragmentCount(int count)
    {
        mFragmentCount += count;
    }

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void addPositionCount(int seIndex) { ++mPositionCounts[seIndex]; }
    public void addPositionCount(int seIndex, int count) { mPositionCounts[seIndex] += count; }

    public String[] getTranscriptNames() { return mTranscriptNames; }
    public String[] getBaseContext() { return mBaseContext; }
    public int[] getNearestExonDistance() { return mNearestExonDistance; }

    public void setGeneId(final String geneId) { mGeneId = geneId; }
    public final String getGeneId() { return mGeneId; }

    public static final String ASJ_TRANS_NONE = "NONE";

    public void calcSummaryData(final GeneReadData gene)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<RegionReadData> regions = (se == SE_START) ? mSjStartRegions : mSjEndRegions;
            mTranscriptNames[se] = regions.isEmpty() ? ASJ_TRANS_NONE : generateTranscriptNames(regions);
            mNearestExonDistance[se] = calcNearestExonBoundary(se, gene);
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

    public int calcNearestExonBoundary(int seIndex, final GeneReadData gene)
    {
        if((seIndex == SE_START && !mSjStartRegions.isEmpty()) || (seIndex == SE_END && !mSjEndRegions.isEmpty()))
            return 0;

        /* find distance to nearest splice acceptor or donor as follows:
            - 5' in intron - go back to last exon end
            - 5' in exon - go forward to exon end - record as -ve value
            - 3' in intron - go forward to next exon start
            - 3' in exon - go back to exon start, record as -ve value
        */

        int position = SpliceJunction[seIndex];
        boolean forwardStrand = gene.GeneData.Strand == 1;
        boolean isFivePrime = (seIndex == SE_START) == forwardStrand;
        boolean isExonic = RegionContexts[seIndex].equals(EXONIC) || RegionContexts[seIndex].equals(MIXED);
        boolean searchForwards = (isFivePrime && isExonic) || (!isFivePrime && !isExonic);

        int nearestBoundary = 0;

        for(RegionReadData region : gene.getExonRegions())
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

    public void setBaseContext(final IndexedFastaSequenceFile refGenome, final String chromosome)
    {
        setJunctionBaseContext(refGenome, new String[] { chromosome, chromosome }, SpliceJunction, mBaseContext);
    }

    public void setDonorAcceptorBases(final String[] baseContext)
    {
        mDonorAcceptorBases = formDonorAcceptorBases(baseContext);
    }

    public String getDonorAcceptorBases() { return mDonorAcceptorBases; }

    public int getKnownSpliceBaseStrand()
    {
        if(mDonorAcceptorBases.equals(SP_SEQ_POS_STRAND_1) || mDonorAcceptorBases.equals(SP_SEQ_POS_STRAND_2))
            return 1;
        else if(mDonorAcceptorBases.equals(SP_SEQ_NEG_STRAND_1) || mDonorAcceptorBases.equals(SP_SEQ_NEG_STRAND_2))
            return -1;
        else
            return 0;
    }

    public void cullNonMatchedTranscripts(final List<Integer> validTransIds)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<RegionReadData> regions = se == SE_START ? mSjStartRegions : mSjEndRegions;

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
        final List<Integer> validTransIds = Lists.newArrayList();

        // if any splice junctions have been matched, restrict the set of transcripts (and by implication genes) to these only
        if(!mSjStartRegions.isEmpty() || !mSjEndRegions.isEmpty())
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                List<RegionReadData> sjRegions = se == SE_START ? mSjStartRegions : mSjEndRegions;

                for (RegionReadData region : sjRegions)
                {
                    region.getTransExonRefs().stream().forEach(x -> validTransIds.add(x.TransId));
                }
            }
        }
        else
        {
            for (RegionReadData region : candidateRegions)
            {
                if (positionWithin(SpliceJunction[SE_START], region.start(), region.end())
                        && positionWithin(SpliceJunction[SE_END], region.start(), region.end()))
                {
                    validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId).collect(Collectors.toList()));
                    continue;
                }

                if (positionWithin(SpliceJunction[SE_START], region.start(), region.end()))
                {
                    // each transcript must be present in the next region to be valid
                    validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId)
                            .filter(x -> region.getPostRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList()));
                }

                if (positionWithin(SpliceJunction[SE_END], region.start(), region.end()))
                {
                    validTransIds.addAll(region.getTransExonRefs().stream().map(x -> x.TransId)
                            .filter(x -> region.getPreRegions().stream().anyMatch(y -> y.hasTransId(x))).collect(Collectors.toList()));
                }
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
                mGeneId != null ? mGeneId : "unset", SpliceJunction[SE_START], SpliceJunction[SE_END],
                RegionContexts[SE_START], RegionContexts[SE_END], mType, mFragmentCount);
    }

    public static String csvHeader()
    {
        return "GeneId,GeneName,Chromosome,Strand,SjStart,SjEnd,FragCount,DepthStart,DepthEnd"
                + ",Type,ContextStart,ContextEnd,NearestStartExon,NearestEndExon"
                + ",BasesStart,BasesEnd,TransStart,TransEnd,OverlappingGenes";
    }

    public String toCsv(final EnsemblGeneData geneData)
    {
        return new StringJoiner(DELIMITER)
                .add(getGeneId())
                .add(geneData.GeneName)
                .add(Chromosome)
                .add(String.valueOf(geneData.Strand))
                .add(String.valueOf(SpliceJunction[SE_START]))
                .add(String.valueOf(SpliceJunction[SE_END]))
                .add(String.valueOf(getFragmentCount()))
                .add(String.valueOf(getPositionCount(SE_START)))
                .add(String.valueOf(getPositionCount(SE_END)))
                .add(String.valueOf(type()))
                .add(String.valueOf(RegionContexts[SE_START]))
                .add(String.valueOf(RegionContexts[SE_END]))
                .add(String.valueOf(getNearestExonDistance()[SE_START]))
                .add(String.valueOf(getNearestExonDistance()[SE_END]))
                .add(getBaseContext()[SE_START])
                .add(getBaseContext()[SE_END])
                .add(getTranscriptNames()[SE_START])
                .add(getTranscriptNames()[SE_END])
                .toString();
    }

    public static AltSpliceJunction fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(DELIMITER);

        final int[] spliceJunction =
                { Integer.parseInt(items[fieldIndexMap.get("SjStart")]), Integer.parseInt(items[fieldIndexMap.get("SjEnd")]) };

        final AltSpliceJunctionContext[] contexts =
                { AltSpliceJunctionContext.valueOf(items[fieldIndexMap.get("ContextStart")]),
                        AltSpliceJunctionContext.valueOf(items[fieldIndexMap.get("ContextEnd")]) };

        AltSpliceJunction altSJ = new AltSpliceJunction(
                items[fieldIndexMap.get("Chromosome")], spliceJunction, AltSpliceJunctionType.valueOf(items[fieldIndexMap.get("Type")]),
                contexts, Lists.newArrayList(), Lists.newArrayList());

        altSJ.setGeneId(items[fieldIndexMap.get("GeneId")]);
        altSJ.addFragmentCount(Integer.parseInt(items[fieldIndexMap.get("FragCount")]));
        altSJ.addPositionCount(SE_START, Integer.parseInt(items[fieldIndexMap.get("DepthStart")]));
        altSJ.addPositionCount(SE_END, Integer.parseInt(items[fieldIndexMap.get("DepthEnd")]));

        final String[] baseContext = { items[fieldIndexMap.get("BasesStart")], items[fieldIndexMap.get("BasesEnd")] };
        altSJ.setDonorAcceptorBases(baseContext);

        altSJ.mTranscriptNames[SE_START] = items[fieldIndexMap.get("TransStart")];
        altSJ.mTranscriptNames[SE_END] = items[fieldIndexMap.get("TransEnd")];

        return altSJ;
    }
}
