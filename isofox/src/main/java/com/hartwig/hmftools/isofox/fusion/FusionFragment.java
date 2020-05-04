package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalAcceptor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalDonor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.JUNCTION_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocationPair;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.lowerChromosome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class FusionFragment
{
    private final List<ReadRecord> mReads;

    private final int[] mGeneCollections;
    private final String[] mChromosomes;
    private final byte[] mOrientations;
    private final int[] mJunctionPositions; // fusion junction if exists
    private final byte[] mJunctionOrientations; // orientation at junction if exists
    private final FusionJunctionType[] mJunctionTypes;
    private final String[] mJunctionBases; // the 10 bases leading up to the junction if it exists
    private boolean mJunctionBasesMatched;
    private final String[] mJunctionSpliceBases; // the 2 donor/acceptor bases
    private FusionFragmentType mType;

    private final RegionMatchType[] mRegionMatchTypes; // top-ranking region match type from the reads
    private final List<TransExonRef>[] mTransExonRefs;

    public FusionFragment(final List<ReadRecord> reads)
    {
        mReads = reads;

        mGeneCollections = new int[SE_PAIR];
        mJunctionPositions = new int[] {-1, -1};
        mChromosomes = new String[] {"", ""};
        mJunctionOrientations = new byte[] {0, 0};
        mOrientations = new byte[] {0, 0};
        mRegionMatchTypes = new RegionMatchType[] { RegionMatchType.NONE, RegionMatchType.NONE };
        mJunctionTypes = new FusionJunctionType[] { FusionJunctionType.UNKNOWN, FusionJunctionType.UNKNOWN };
        mJunctionBases = new String[] {"", ""};
        mJunctionBasesMatched = false;
        mJunctionSpliceBases = new String[] {"", ""};

        mTransExonRefs = new List[SE_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();

        // divide reads into the 2 gene collections
        final List<String> chrGeneCollections = Lists.newArrayListWithCapacity(2);
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);
        final List<Integer> positions = Lists.newArrayListWithCapacity(2);
        final List<Byte> orientations = Lists.newArrayListWithCapacity(2);
        final Map<String,List<ReadRecord>> readGroups = Maps.newHashMap();

        for(final ReadRecord read : reads)
        {
            final String chrGeneId = formLocation(read.Chromosome, read.getGeneCollecton());

            List<ReadRecord> readGroup = readGroups.get(chrGeneId);

            if(readGroup == null)
            {
                readGroups.put(chrGeneId, Lists.newArrayList(read));

                chrGeneCollections.add(chrGeneId);
                chromosomes.add(read.Chromosome);
                positions.add(read.PosStart); // no overlap in gene collections so doesn't matter which position is used
                orientations.add(read.orientation());
            }
            else
            {
                readGroup.add(read);
            }
        }

        boolean overlappingReads = hasOverlappingReadGroups(readGroups);

        if(overlappingReads || readGroups.size() > 2)
        {
            mChromosomes[SE_START] = mChromosomes[SE_END] = chromosomes.get(0);

            mGeneCollections[SE_START] = reads.stream().mapToInt(x -> x.getGeneCollecton()).min().orElse(-1);
            mGeneCollections[SE_END] = reads.stream().mapToInt(x -> x.getGeneCollecton()).max().orElse(-1);
            mType = UNKNOWN;
            return;
        }

        // first determine which is the start and end chromosome & position as for SVs
        if(readGroups.size() > 1)
        {
            int lowerIndex;

            if (chromosomes.get(0).equals(chromosomes.get(1)))
                lowerIndex = positions.get(0) <= positions.get(1) ? 0 : 1;
            else
                lowerIndex = lowerChromosome(chromosomes.get(0), chromosomes.get(1)) ? 0 : 1;

            for (int se = SE_START; se <= SE_END; ++se)
            {
                int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
                final String chrGeneId = chrGeneCollections.get(index);

                final List<ReadRecord> readGroup = readGroups.get(chrGeneId);

                mChromosomes[se] = chromosomes.get(index);
                mGeneCollections[se] = readGroup.get(0).getGeneCollecton();

                if(readGroup.size() == 2)
                    mOrientations[se] = readGroup.stream()
                            .filter(x -> x.getSuppAlignment() == null).findFirst().map(x -> x.orientation()).orElse((byte)0);
                else
                    mOrientations[se] = readGroup.get(0).orientation();
            }

            if(mReads.stream().anyMatch(x -> x.getSuppAlignment() != null))
                setJunctionData(readGroups, chrGeneCollections, lowerIndex);

            extractTranscriptExonData(true);

            mType = calcType();
        }
        else
        {
            // these fragments are from secondary searches for one-side support of the fusion, and so will only have one gene collection
            // and will be missing transcript and exon information
            mChromosomes[SE_START] = mChromosomes[SE_END] = chromosomes.get(0);
            mGeneCollections[SE_START] = mGeneCollections[SE_END] = reads.get(0).getGeneCollecton();
            mOrientations[SE_START] = mOrientations[SE_END] = reads.get(0).orientation();

            setJunctionData(readGroups, chrGeneCollections, 0);
            extractTranscriptExonData(false);

            mType = UNKNOWN;
        }
    }

    public String readId() { return mReads.get(0).Id; }

    public final List<ReadRecord> getReads() { return mReads; }
    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final byte[] orientations() { return mOrientations; }

    public boolean isSingleGene()
    {
        return mChromosomes[SE_START].equals(mChromosomes[SE_END]) && mGeneCollections[SE_START] == mGeneCollections[SE_END];
    }

    public final int[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }
    public final String[] junctionBases() { return mJunctionBases; }

    public final RegionMatchType[] regionMatchTypes() { return mRegionMatchTypes; }
    public final FusionJunctionType[] junctionTypes() { return mJunctionTypes; }

    public boolean isUnspliced() { return mRegionMatchTypes[SE_START] == INTRON && mRegionMatchTypes[SE_END] == INTRON; }
    public boolean isSpliced() { return exonBoundary(mRegionMatchTypes[SE_START]) && exonBoundary(mRegionMatchTypes[SE_END]); }

    public String locationPair() { return formLocationPair(mChromosomes, mGeneCollections); }

    public static boolean isRealignedFragmentCandidate(final ReadRecord read)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if (read.isSoftClipped(se))
            {
                int scLength =
                        se == SE_START ? read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

                if (scLength >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH && scLength <= REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                    return true;
            }
        }

        return false;
    }

    private boolean hasOverlappingReadGroups(final Map<String,List<ReadRecord>> readGroups)
    {
        if(readGroups.size() == 1)
            return false;

        for(Map.Entry<String,List<ReadRecord>> entry1 : readGroups.entrySet())
        {
            for(Map.Entry<String,List<ReadRecord>> entry2 : readGroups.entrySet())
            {
                if(entry1.getKey().equals(entry2.getKey()))
                    continue;

                if(entry1.getValue().stream().anyMatch(x -> entry2.getValue().stream()
                        .anyMatch(y -> positionsOverlap(x.PosStart, x.PosEnd, y.PosStart, y.PosEnd))))
                {
                    return true;
                }
            }
        }

        return false;
    }

    private void setJunctionData(final Map<String,List<ReadRecord>> readGroups, final List<String> chrGeneCollections, int lowerIndex)
    {
        // find the reads with supplementary read info and use this to set
        final String[] softClipBases = new String[] {"", ""};
        boolean[] junctionValid = { false, false };

        boolean isSingleGene = readGroups.size() == 1;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(isSingleGene && se == SE_END)
                return;

            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);

            // find the outermost soft-clipped read to use for the splice junction position
            int sjPosition = 0;
            byte sjOrientation = 0;

            final List<ReadRecord> readGroup = readGroups.get(chrGeneId);

            ReadRecord read = null;
            int maxScLength = 0;

            for(ReadRecord rgRead : readGroup)
            {
                if(!rgRead.Cigar.containsOperator(CigarOperator.S))
                    continue;

                if(!isSingleGene)
                {
                    if(rgRead.getSuppAlignment() != null)
                        read = rgRead;
                }
                else
                {
                    int scLeft = rgRead.isSoftClipped(SE_START) ? rgRead.Cigar.getFirstCigarElement().getLength() : 0;
                    int scRight = rgRead.isSoftClipped(SE_END) ? rgRead.Cigar.getLastCigarElement().getLength() : 0;

                    if(max(scLeft, scRight) > maxScLength)
                    {
                        maxScLength = max(scLeft, scRight);
                        read = rgRead;
                    }
                }
            }

            if(read == null)
                continue;

            int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
            int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

            boolean useLeft = false;

            if(scLeft > 0 && scRight > 0)
            {
                // should be very unlikely since implies a very short exon and even then would expect it to be mapped
                if(scLeft >= scRight)
                {
                    useLeft = true;
                }
                else if(scRight > scLeft)
                {
                    useLeft = false;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                useLeft = scLeft > 0;
            }

            if(useLeft)
            {
                sjPosition = read.getCoordsBoundary(SE_START);
                sjOrientation = -1;
            }
            else
            {
                sjPosition = read.getCoordsBoundary(SE_END);
                sjOrientation = 1;
            }

            mJunctionPositions[se] = sjPosition;
            mJunctionOrientations[se] = sjOrientation;
            junctionValid[se] = true;

            /*  Junction bases will now be set from the ref genome, not the read, in case they contain SNVs and cause invalid realignments

            int baseLength = min(10, read.Length);

            if(mJunctionOrientations[se] == 1)
            {
                int scLength = read.Cigar.getLastCigarElement().getLength();
                int readEndPos = read.Length - scLength;
                mJunctionBases[se] = read.ReadBases.substring(readEndPos - baseLength, readEndPos);
                // softClipBases[se] = read.ReadBases.substring(readEndPos, readEndPos + baseLength);
            }
            else
            {
                int scLength = read.Cigar.getFirstCigarElement().getLength();
                int readStartPos = scLength;
                mJunctionBases[se] = read.ReadBases.substring(readStartPos, readStartPos + baseLength);
                // softClipBases[se] = read.ReadBases.substring(readStartPos - baseLength, readStartPos);
            }
            */
        }

        if(junctionValid[SE_START] && junctionValid[SE_END])
            mJunctionBasesMatched = true;

        /*
        if(!mJunctionBases[SE_START].isEmpty() && !mJunctionBases[SE_END].isEmpty())
        {
            if(mJunctionBases[SE_START].equals(softClipBases[SE_END]) && mJunctionBases[SE_END].equals(softClipBases[SE_START]))
            {
                mJunctionBasesMatched = true;
            }
        }

\       */
    }

    public void setType(FusionFragmentType type) { mType = type; }

    private FusionFragmentType calcType()
    {
        if(mJunctionBasesMatched)
        {
            return MATCHED_JUNCTION;
        }
        else
        {
            return DISCORDANT;
        }
    }

    public void setGeneData(int geneCollection, final List<TransExonRef> transExonRefs)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            mGeneCollections[se] = geneCollection;
            mTransExonRefs[se].addAll(transExonRefs);
        }
    }

    public StructuralVariantType getImpliedSvType()
    {
        return impliedSvType(mChromosomes, mJunctionOrientations);
    }

    public final List<TransExonRef>[] getTransExonRefs() { return mTransExonRefs; }

    public List<String> getGeneIds(int seIndex)
    {
        final List<String> geneIds = Lists.newArrayList();

        for(TransExonRef transExonRef : mTransExonRefs[seIndex])
        {
            if(!geneIds.contains(transExonRef.GeneId))
                geneIds.add(transExonRef.GeneId);
        }

        return geneIds;
    }

    private void extractTranscriptExonData(boolean checkJunction)
    {
        // take the transcript & exon refs from each read
        for(int se = SE_START; se <= SE_END; ++se)
        {
            for (final ReadRecord read : mReads)
            {
                if (!read.Chromosome.equals(mChromosomes[se]) || read.getGeneCollecton() != mGeneCollections[se])
                    continue;

                if(checkJunction && mJunctionPositions[se] > 0 && !positionWithin(mJunctionPositions[se], read.PosStart, read.PosEnd))
                    continue;

                for (Map.Entry<RegionMatchType, List<TransExonRef>> entry : read.getTransExonRefs().entrySet())
                {
                    RegionMatchType matchType = entry.getKey();

                    if(matchRank(matchType) > matchRank(mRegionMatchTypes[se]))
                        mRegionMatchTypes[se] = matchType;

                    for(TransExonRef readTransExonRef : entry.getValue())
                    {
                        boolean found = false;

                        for(TransExonRef transExonRef : mTransExonRefs[se])
                        {
                            if (transExonRef.TransId == readTransExonRef.TransId)
                            {
                                found = true;

                                if (transExonRef.ExonRank != readTransExonRef.ExonRank)
                                {
                                    ISF_LOGGER.trace("multi-exon: read({} cigar={}) ref1({}) ref2({})",
                                            read.Id, read.Cigar.toString(), transExonRef.toString(), readTransExonRef.toString());

                                    // will be handled later on
                                    mTransExonRefs[se].add(readTransExonRef);
                                }

                                break;
                            }
                        }

                        if(!found)
                        {
                            mTransExonRefs[se].add(readTransExonRef);
                        }
                    }
                }
            }
        }
    }

    public void validateTranscriptExons(final List<TranscriptData> transDataList, int seIndex)
    {
        if(mJunctionPositions[seIndex] < 0 || !exonBoundary(mRegionMatchTypes[seIndex]))
            return;

        int junctionPosition = mJunctionPositions[seIndex];

        int index = 0;
        while(index < mTransExonRefs[seIndex].size())
        {
            TransExonRef transExonRef = mTransExonRefs[seIndex].get(index);
            TranscriptData transData = transDataList.stream().filter(x -> x.TransId == transExonRef.TransId).findFirst().orElse(null);

            if(transData == null)
            {
                mTransExonRefs[seIndex].remove(index);
                continue;
            }

            ExonData exon = transData.exons().stream()
                    .filter(x -> x.ExonStart == junctionPosition || x.ExonEnd == junctionPosition).findFirst().orElse(null);

            if(exon != null && transExonRef.ExonRank != exon.ExonRank)
            {
                mTransExonRefs[seIndex].remove(index);
                continue;
            }

            ++index;
        }
    }

    public void setJunctionTypes(final IndexedFastaSequenceFile refGenome, final byte[] junctionStrands)
    {
        if(refGenome == null)
            return;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if (exonBoundary(mRegionMatchTypes[se]))
            {
                mJunctionTypes[se] = FusionJunctionType.KNOWN;
            }
            else if(mJunctionPositions[se] > 0)
            {
                String daBases = mJunctionSpliceBases[se];

                if(junctionStrands != null)
                {
                    boolean isDonor = (mJunctionOrientations[se] == junctionStrands[se]);

                    if (isDonor && canonicalDonor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                    else if (!isDonor && canonicalAcceptor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                }
                else
                {
                    // try them both
                    byte asDonorStrand = mJunctionOrientations[se];
                    byte asAcceptorStrand = (byte)(-mJunctionOrientations[se]);

                    if (canonicalDonor(daBases, asDonorStrand) || canonicalAcceptor(daBases, asAcceptorStrand))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                }
            }
        }
    }

    public void setJunctionBases(final IndexedFastaSequenceFile refGenome)
    {
        if(mType != MATCHED_JUNCTION)
            return;

        if(refGenome != null)
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                int junctionBase = mJunctionPositions[se];

                if (junctionOrientations()[se] == 1)
                {
                    String junctionBases = refGenome.getSubsequenceAt(
                            mChromosomes[se], junctionBase - JUNCTION_BASE_LENGTH + 1, junctionBase + 2).getBaseString();

                    mJunctionBases[se] = junctionBases.substring(0, JUNCTION_BASE_LENGTH);
                    mJunctionSpliceBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH);
                }
                else
                {
                    String junctionBases = refGenome.getSubsequenceAt(
                            mChromosomes[se], junctionBase - 2, junctionBase + JUNCTION_BASE_LENGTH - 1).getBaseString();

                    mJunctionBases[se] = junctionBases.substring(2);
                    mJunctionSpliceBases[se] = junctionBases.substring(0, 2);
                }
            }
        }
        else
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                final int seIndex = se;
                int junctionBase = mJunctionPositions[se];

                ReadRecord read = mReads.stream()
                        .filter(x -> x.Chromosome.equals(mChromosomes[seIndex]))
                        .filter(x -> x.getCoordsBoundary(junctionOrientations()[seIndex] == 1 ? SE_END : SE_START) == junctionBase)
                        .findFirst().orElse(null);

                if (junctionOrientations()[se] == 1)
                {
                    mJunctionBases[se] = read.ReadBases.substring(read.Length - JUNCTION_BASE_LENGTH, read.Length);
                }
                else
                {
                    mJunctionBases[se] = read.ReadBases.substring(0, JUNCTION_BASE_LENGTH);
                }
            }
        }
    }

}
