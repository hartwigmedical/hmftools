package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_ACCEPTOR;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_DONOR_1;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_DONOR_2;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_NEG_STRAND_ACCEPTOR;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_NEG_STRAND_DONOR_1;
import static com.hartwig.hmftools.isofox.common.CommonUtils.SP_SEQ_NEG_STRAND_DONOR_2;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.EXONIC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.MIXED;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.CIRCULAR;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

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
    private String mGeneName;

    private final String[] mSelectedTranscripts;
    private final int[] mSelectedExons;

    // calculated values
    private final String[] mTranscriptNames;
    private final String[] mBaseContext;
    private String mDonorAcceptorBases;
    private final int[] mNearestExonDistance;

    private int mCohortFrequency;
    private final String mInitialReadId;
    private String mFilter;

    public AltSpliceJunction(
            final String chromosome, final int[] spliceJunction, AltSpliceJunctionType type, final String initialReadId,
            final AltSpliceJunctionContext[] regionContexts, final List<RegionReadData> sjStartRegions, final List<RegionReadData> sjEndRegions)
    {
        Chromosome = chromosome;
        SpliceJunction = spliceJunction;
        RegionContexts = regionContexts;

        mSjStartRegions = sjStartRegions;
        mSjEndRegions = sjEndRegions;

        mCandidateTransIds = Lists.newArrayList();
        mGeneId = null;
        mGeneName = null;
        mInitialReadId = initialReadId;

        mType = type;

        mFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
        mTranscriptNames = new String[SE_PAIR];
        mBaseContext = new String[SE_PAIR];
        mNearestExonDistance = new int[SE_PAIR];
        mDonorAcceptorBases = "";

        mSelectedTranscripts = new String[SE_PAIR];
        mSelectedExons = new int[] {-1, -1};
        mCohortFrequency = -1;
        mFilter = "";
    }

    public boolean matches(final AltSpliceJunction other)
    {
        return Chromosome.equals(other.Chromosome)
                && SpliceJunction[SE_START] == other.SpliceJunction[SE_START]
                && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public AltSpliceJunctionType type() { return mType; }
    public void overrideType(AltSpliceJunctionType type) { mType = type; }

    public int length() { return SpliceJunction[SE_END] - SpliceJunction[SE_START]; }

    public final List<RegionReadData> getSjStartRegions() { return mSjStartRegions; }
    public final List<RegionReadData> getSjEndRegions() { return mSjEndRegions; }

    public int getFragmentCount() { return mFragmentCount;}
    public void addFragmentCount() { ++mFragmentCount;}

    public void addPositionCount(int seIndex, int count) { mPositionCounts[seIndex] += count; }

    public void setGeneData(final String geneId, final String geneName)
    {
        mGeneId = geneId;
        mGeneName = geneName;
    }

    public String geneId() { return mGeneId; }
    public String geneName() { return mGeneName; }

    public int cohortFrequency() { return mCohortFrequency; }
    public void setCohortFrequency(int frequency) { mCohortFrequency = frequency; }

    public String filter() { return mFilter; }
    public void setFilter(final String filter) { mFilter = filter; }

    public static final String ASJ_TRANS_NONE = "NONE";

    public void calcSummaryData(final GeneReadData gene)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<RegionReadData> regions = (se == SE_START) ? mSjStartRegions : mSjEndRegions;
            mTranscriptNames[se] = regions.isEmpty() ? ASJ_TRANS_NONE : generateTranscriptNames(regions);
            mNearestExonDistance[se] = calcNearestExonBoundary(se, gene);
        }
    }

    private String generateTranscriptNames(final List<RegionReadData> regions)
    {
        List<Integer> validTransIds = candidateTransIds();
        StringJoiner transNames = new StringJoiner(";");

        for(RegionReadData region: regions)
        {
            region.getTransExonRefs().stream()
                    .filter(x -> validTransIds.contains(x.TransId))
                    .forEach(x -> transNames.add(x.TransName));
        }

        return transNames.toString();
    }

    private void setTranscriptAndExon(final String[] transcripts, final int[] exons, final int seIndex, final List<RegionReadData> regions)
    {
        if(regions.isEmpty())
            return;

        for(RegionReadData regionReadData : regions)
        {
            for(TransExonRef transExonRef : regionReadData.getTransExonRefs())
            {
                if(transExonRef.isCanonical())
                {
                    transcripts[seIndex] = transExonRef.TransName;
                    exons[seIndex] = transExonRef.ExonRank;
                    return;
                }
            }
        }

        // set to any - could find one with most support
        if(!regions.get(0).getTransExonRefs().isEmpty())
        {
            TransExonRef transExonRef = regions.get(0).getTransExonRefs().get(0);
            transcripts[seIndex] = transExonRef.TransName;
            exons[seIndex] = transExonRef.ExonRank;
        }
    }

    public int calcNearestExonBoundary(int seIndex, final GeneReadData gene)
    {
        if((seIndex == SE_START && !mSjStartRegions.isEmpty()) || (seIndex == SE_END && !mSjEndRegions.isEmpty()))
        {
            // an exact splice junction match - set transcript and exon, favouring canonical if matched
            setTranscriptAndExon(mSelectedTranscripts, mSelectedExons, seIndex, seIndex == SE_START ? mSjStartRegions : mSjEndRegions);
            return 0;
        }

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
        RegionReadData nearestRegion = null;

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
                distance = searchForwards ? position - region.end() : region.start() - position;

                if(nearestBoundary == 0 || distance > nearestBoundary)
                {
                    nearestBoundary = distance;
                    nearestRegion = region;
                }
            }
            else
            {
                if(isExonic)
                    continue;

                if((searchForwards && region.start() < position) || (!searchForwards && position < region.end()))
                    continue;

                // will be positive
                distance = searchForwards ? region.start() - position : position - region.end();

                if(nearestBoundary == 0 || distance < nearestBoundary)
                {
                    nearestBoundary = distance;
                    nearestRegion = region;
                }
            }
        }

        if(nearestRegion != null)
        {
            setTranscriptAndExon(mSelectedTranscripts, mSelectedExons, seIndex, List.of(nearestRegion));
        }

        return nearestBoundary;
    }

    public void setBaseContext(final RefGenomeInterface refGenome, final String chromosome)
    {
        // get the 2 bases leading up to and including the splice junction base, and 10 bases into the unspliced region
        // the donor/acceptor bases are the 2 bases leading up to the junction from the unspliced side
        for(int se = SE_START; se <= SE_END; ++se)
        {
            boolean startInExon = mType != CIRCULAR; // DUP-like SJs are reversed

            int startOffset = ((se == SE_START) == startInExon) ? 1 : 10;
            int endOffset = startOffset == 1 ? 10: 1;

            mBaseContext[se] = refGenome.getBaseString(
                    chromosome, SpliceJunction[se] - startOffset, SpliceJunction[se] + endOffset);
        }
    }

    private static final String DA_DELIM = "-";

    public static final String SP_SEQ_POS_STRAND_1 = SP_SEQ_DONOR_1 + DA_DELIM + SP_SEQ_ACCEPTOR;// "GT-AG";
    public static final String SP_SEQ_POS_STRAND_2 = SP_SEQ_DONOR_2 + DA_DELIM + SP_SEQ_ACCEPTOR; // "GC-AG";
    public static final String SP_SEQ_NEG_STRAND_1 = SP_SEQ_NEG_STRAND_ACCEPTOR + DA_DELIM + SP_SEQ_NEG_STRAND_DONOR_1; // "CT-AC";
    public static final String SP_SEQ_NEG_STRAND_2 = SP_SEQ_NEG_STRAND_ACCEPTOR + DA_DELIM + SP_SEQ_NEG_STRAND_DONOR_2; // "CT-GC";

    private static String startDonorAcceptorBases(final String baseContext)
    {
        return baseContext.length() == 12 ? baseContext.substring(2,4) : "";
    }

    private static String endDonorAcceptorBases(final String baseContext)
    {
        return baseContext.length() == 12 ? baseContext.substring(8,10) : "";
    }

    public static String getDonorAcceptorBases(final String[] baseContext)
    {
        return startDonorAcceptorBases(baseContext[SE_START]) + DA_DELIM + endDonorAcceptorBases(baseContext[SE_END]);
    }

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
                ++index;
            else
                mCandidateTransIds.remove(index);
        }
    }

    public List<Integer> candidateTransIds() { return mCandidateTransIds; }

    public void setCandidateTranscripts(final List<RegionReadData> candidateRegions)
    {
        List<Integer> validTransIds = Lists.newArrayList();

        // if any splice junctions have been matched, restrict the set of transcripts (and by implication genes) to these only
        if(!mSjStartRegions.isEmpty() || !mSjEndRegions.isEmpty())
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                List<RegionReadData> sjRegions = se == SE_START ? mSjStartRegions : mSjEndRegions;

                for(RegionReadData region : sjRegions)
                {
                    region.getTransExonRefs().stream().forEach(x -> validTransIds.add(x.TransId));
                }
            }
        }
        else
        {
            for(RegionReadData region : candidateRegions)
            {
                if(positionWithin(SpliceJunction[SE_START], region.start(), region.end())
                        && positionWithin(SpliceJunction[SE_END], region.start(), region.end()))
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

    // file output and conversion
    public NovelSpliceJunction convert()
    {
        return ImmutableNovelSpliceJunction.builder()
                .geneName(mGeneName)
                .chromosome(Chromosome)
                .junctionStart(SpliceJunction[SE_START])
                .junctionEnd(SpliceJunction[SE_END])
                .type(mType)
                .transcriptStart(mSelectedTranscripts[SE_START] != null ? mSelectedTranscripts[SE_START] : "")
                .transcriptEnd(mSelectedTranscripts[SE_END] != null ? mSelectedTranscripts[SE_END] : "")
                .exonStart(mSelectedExons[SE_START])
                .exonEnd(mSelectedExons[SE_END])
                .fragmentCount(mFragmentCount)
                .depthStart(mPositionCounts[SE_START])
                .depthEnd(mPositionCounts[SE_END])
                .regionStart(RegionContexts[SE_START])
                .regionEnd(RegionContexts[SE_END])
                .basesStart(mBaseContext[SE_START])
                .basesEnd(mBaseContext[SE_END])
                .cohortFrequency(mCohortFrequency)
                .build();
    }

    public static final String FLD_NEAREST_EXON_START = "NearestStartExon";
    public static final String FLD_NEAREST_EXON_END = "NearestEndExon";
    public static final String FLD_INIT_READ_ID = "InitialReadId";
    public static final String FLD_TRANS_START = "TransStart";
    public static final String FLD_TRANS_END = "TransEnd";
    public static final String FLD_FILTER = "Filter";

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(FLD_GENE_ID);
        sj.add(NovelSpliceJunctionFile.header());

        sj.add(FLD_FILTER);
        sj.add(FLD_TRANS_START);
        sj.add(FLD_TRANS_END);
        sj.add(FLD_NEAREST_EXON_START);
        sj.add(FLD_NEAREST_EXON_END);
        sj.add(FLD_INIT_READ_ID);

        return sj.toString();
    }

    public String write(final NovelSpliceJunction novelSpliceJunction)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(mGeneId);
        sj.add(NovelSpliceJunctionFile.write(novelSpliceJunction));
        sj.add(mFilter.isEmpty() ? "PASS" : mFilter);
        sj.add(mTranscriptNames[SE_START]);
        sj.add(mTranscriptNames[SE_END]);
        sj.add(String.valueOf(mNearestExonDistance[SE_START]));
        sj.add(String.valueOf(mNearestExonDistance[SE_END]));
        sj.add(mInitialReadId);

        return sj.toString();
    }
}
