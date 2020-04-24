package com.hartwig.hmftools.isofox.fusion;

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
import static com.hartwig.hmftools.isofox.common.RnaUtils.endDonorAcceptorBases;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.startDonorAcceptorBases;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.ONE_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.formLocationPair;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.lowerChromosome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RnaUtils;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class FusionFragment
{
    private final List<ReadRecord> mReads;

    private final int[] mGeneCollections;
    private final String[] mChromosomes;
    private final long[] mJunctionPositions; // fusion junction is exists
    private final byte[] mJunctionOrientations; // orientation at junction
    private final boolean[] mJunctionValid;
    private final FusionJunctionType[] mJunctionTypes;
    private final String[] mJunctionBaseContext;
    private FusionFragmentType mType;

    private final RegionMatchType[] mRegionMatchTypes; // top-ranking region match type from the reads
    private final List<TransExonRef>[] mTransExonRefs;

    public FusionFragment(final List<ReadRecord> reads)
    {
        mReads = reads;

        mGeneCollections = new int[SE_PAIR];
        mJunctionPositions = new long[] {-1, -1};
        mChromosomes = new String[] {"", ""};
        mJunctionOrientations = new byte[] {0, 0};
        mJunctionValid = new boolean[] {false, false};
        mRegionMatchTypes = new RegionMatchType[] { RegionMatchType.NONE, RegionMatchType.NONE };
        mJunctionTypes = new FusionJunctionType[] { FusionJunctionType.UNKNOWN, FusionJunctionType.UNKNOWN };
        mJunctionBaseContext = new String[] {"", ""};

        mTransExonRefs = new List[SE_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();

        // divide reads into the 2 gene collections
        final List<String> chrGeneCollections = Lists.newArrayListWithCapacity(2);
        final List<String> chromosomes = Lists.newArrayListWithCapacity(2);
        final List<Long> positions = Lists.newArrayListWithCapacity(2);
        final Map<String,List<ReadRecord>> readGroups = Maps.newHashMap();

        for(final ReadRecord read : reads)
        {
            final String chrGeneId = read.chromosomeGeneId();

            List<ReadRecord> readGroup = readGroups.get(chrGeneId);

            if(readGroup == null)
            {
                readGroups.put(chrGeneId, Lists.newArrayList(read));

                chrGeneCollections.add(chrGeneId);
                chromosomes.add(read.Chromosome);
                positions.add(read.PosStart); // no overlap in gene collections so doesn't matter which position is used
            }
            else
            {
                readGroup.add(read);
            }
        }

        // first determine which is the start and end chromosome & position as for SVs
        int lowerIndex;

        if(chromosomes.get(0).equals(chromosomes.get(1)))
            lowerIndex = positions.get(0) < positions.get(1) ? 0 : 1;
        else
            lowerIndex = lowerChromosome(chromosomes.get(0), chromosomes.get(1)) ? 0 : 1;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int index = se == SE_START ? lowerIndex : switchIndex(lowerIndex);
            final String chrGeneId = chrGeneCollections.get(index);

            // find the outermost soft-clipped read to use for the splice junction position
            long sjPosition = 0;
            byte sjOrientation = 0;
            int maxSoftClipping = 0;

            final List<ReadRecord> readGroup = readGroups.get(chrGeneId);

            mChromosomes[se] = chromosomes.get(index);
            mGeneCollections[se] = readGroup.get(0).getGeneCollecton();

            for(ReadRecord read : readGroup)
            {
                if(!read.Cigar.containsOperator(CigarOperator.S))
                    continue;

                int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
                int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

                boolean useLeft = false;

                if(scLeft > 0 && scRight > 0)
                {
                    // should be very unlikely since implies a very short exon and even then would expect it to be mapped
                    if(scLeft >= scRight && scLeft > maxSoftClipping)
                    {
                        maxSoftClipping = scLeft;
                        useLeft = true;
                    }
                    else if(scRight > scLeft && scRight > maxSoftClipping)
                    {
                        maxSoftClipping = scRight;
                        useLeft = false;
                    }
                    else
                    {
                        continue;
                    }
                }
                else if(scLeft > maxSoftClipping)
                {
                    maxSoftClipping = scLeft;
                    useLeft = true;
                }
                else if(scRight > maxSoftClipping)
                {
                    maxSoftClipping = scRight;
                    useLeft = false;
                }
                else
                {
                    continue;
                }

                if(useLeft)
                {
                    sjPosition = read.getCoordsBoundary(true);
                    sjOrientation = -1;
                }
                else
                {
                    sjPosition = read.getCoordsBoundary(false);
                    sjOrientation = 1;
                }
            }

            if(maxSoftClipping > 0)
            {
                mJunctionPositions[se] = sjPosition;
                mJunctionOrientations[se] = sjOrientation;
                mJunctionValid[se] = true;
            }
        }

        extractTranscriptExonData();

        mType = calcType();
    }

    public void setType(FusionFragmentType type) { mType = type; }

    private FusionFragmentType calcType()
    {
        if(mJunctionValid[SE_START] && mJunctionValid[SE_END])
        {
            return BOTH_JUNCTIONS;
        }
        else if(mJunctionValid[SE_START] || mJunctionValid[SE_END])
        {
            return ONE_JUNCTION;
        }
        else
        {
            return DISCORDANT;
        }
    }

    public String readId() { return mReads.get(0).Id; }

    public final List<ReadRecord> getReads() { return mReads; }
    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }

    public final long[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }
    public final boolean[] junctionValid() { return mJunctionValid; }
    public boolean hasBothJunctions() { return mJunctionValid[SE_START] && mJunctionValid[SE_END]; }
    public final RegionMatchType[] regionMatchTypes() { return mRegionMatchTypes; }
    public final FusionJunctionType[] junctionTypes() { return mJunctionTypes; }

    public boolean isUnspliced() { return mRegionMatchTypes[SE_START] == INTRON && mRegionMatchTypes[SE_END] == INTRON; }
    public boolean isSpliced() { return exonBoundary(mRegionMatchTypes[SE_START]) && exonBoundary(mRegionMatchTypes[SE_END]); }

    public String locationPair() { return formLocationPair(mChromosomes, mGeneCollections); }

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

    private void extractTranscriptExonData()
    {
        // take the transcript & exon refs from each read
        for(int se = SE_START; se <= SE_END; ++se)
        {
            for (final ReadRecord read : mReads)
            {
                if (!read.Chromosome.equals(mChromosomes[se]) || read.getGeneCollecton() != mGeneCollections[se])
                    continue;

                if(mJunctionValid[se] && !positionWithin(mJunctionPositions[se], read.PosStart, read.PosEnd))
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
        if(!mJunctionValid[seIndex] || !exonBoundary(mRegionMatchTypes[seIndex]))
            return;

        long junctionPosition = mJunctionPositions[seIndex];

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

        if(hasBothJunctions())
        {
            RnaUtils.setJunctionBaseContext(refGenome, mChromosomes, mJunctionPositions, mJunctionBaseContext);
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if (exonBoundary(mRegionMatchTypes[se]))
            {
                mJunctionTypes[se] = FusionJunctionType.KNOWN;
            }
            else if(mJunctionValid[se])
            {
                String daBases = se == SE_START ? startDonorAcceptorBases(mJunctionBaseContext[se]) : endDonorAcceptorBases(mJunctionBaseContext[se]);

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

}
