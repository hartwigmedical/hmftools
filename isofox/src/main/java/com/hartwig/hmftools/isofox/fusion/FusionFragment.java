package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
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
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.CigarOperator;

public class FusionFragment
{
    private final List<ReadRecord> mReads;

    private final int[] mGeneCollections;
    private final String[] mChromosomes;
    private final long[] mSjPositions;
    private final byte[] mSjOrientations;
    private final boolean[] mSjValid;
    private FusionFragmentType mType;

    private final RegionMatchType[] mRegionMatchTypes; // top-ranking region match type from the reads
    private final List<TransExonRef>[] mTransExonRefs;

    public FusionFragment(final List<ReadRecord> reads)
    {
        mReads = reads;

        mGeneCollections = new int[SE_PAIR];
        mSjPositions = new long[] {-1, -1};
        mChromosomes = new String[] {"", ""};
        mSjOrientations = new byte[] {0, 0};
        mSjValid = new boolean[] {false, false};
        mRegionMatchTypes = new RegionMatchType[] { RegionMatchType.NONE, RegionMatchType.NONE };

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

                int scLeft = read.Cigar.isLeftClipped() ? read.Cigar.getFirstCigarElement().getLength() : 0;
                int scRight = read.Cigar.isRightClipped() ? read.Cigar.getLastCigarElement().getLength() : 0;

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
                mSjPositions[se] = sjPosition;
                mSjOrientations[se] = sjOrientation;
                mSjValid[se] = true;
            }
        }

        mType = setType();

        extractTranscriptExonData();
    }

    private FusionFragmentType setType()
    {
        if(mSjValid[SE_START] && mSjValid[SE_END])
        {
            return BOTH_JUNCTIONS;
        }
        else if(mSjValid[SE_START] || mSjValid[SE_END])
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

    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }
    public boolean hasValidSpliceData() { return mSjValid[SE_START] && mSjValid[SE_END]; }
    public final RegionMatchType[] regionMatchTypes() { return mRegionMatchTypes; }

    public boolean isUnspliced() { return mRegionMatchTypes[SE_START] == INTRON && mRegionMatchTypes[SE_END] == INTRON; }
    public boolean isSpliced() { return exonBoundary(mRegionMatchTypes[SE_START]) && exonBoundary(mRegionMatchTypes[SE_END]); }

    public String locationPair() { return formLocationPair(mChromosomes, mGeneCollections); }

    public static boolean validPositions(final long[] position) { return position[SE_START] > 0 && position[SE_END] > 0; }

    public StructuralVariantType getImpliedSvType()
    {
        return impliedSvType(mChromosomes, mSjOrientations);
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

                if(mSjValid[se] && !positionWithin(mSjPositions[se], read.PosStart, read.PosEnd))
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
                                    ISF_LOGGER.warn("FIXME: need to decide on correct exon: read({}) ref1({}) ref2({})",
                                            read.toString(), transExonRef.toString(), readTransExonRef.toString());

                                    // need to know strand to choose the exon closer to the junction if the junction is valid
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

    /*

    public void setSplicedTransExonRefs(int seIndex)
    {
        // each fragment supporting the splice junction will have the same set of candidate genes
        for(final ReadRecord read : mReads)
        {
            if(!read.Chromosome.equals(mChromosomes[seIndex]))
                continue;

            for(Map.Entry<RegionMatchType,List<TransExonRef>> entry : read.getTransExonRefs().entrySet())
            {
                if(entry.getKey() != EXON_BOUNDARY && entry.getKey() != EXON_MATCH)
                    continue;

                if(read.getCoordsBoundary(true) == mSjPositions[seIndex] || read.getCoordsBoundary(false) == mSjPositions[seIndex])
                {
                    for(TransExonRef transExonRef : entry.getValue())
                    {
                        List<TransExonRef> transExonRefs = mTransExonRefs.get(seIndex).get(entry.getKey());

                        if(transExonRefs == null)
                        {
                            mTransExonRefs.get(seIndex).put(entry.getKey(), Lists.newArrayList(transExonRef));
                        }
                        else
                        {
                            transExonRefs.add(transExonRef);
                        }

                    }
                }
            }
        }
    }

    public void populateUnsplicedTransExonRefs(final List<TranscriptData> transDataList, int seIndex)
    {
        if(!mSjValid[seIndex])
            return;

        long position = mSjPositions[seIndex];
        byte orientation = mSjOrientations[seIndex];

        for (ReadRecord read : mReads)
        {
            if (!read.Chromosome.equals(mChromosomes[seIndex]))
                continue;

            final List<TransExonRef> transExonRefs = read.getTransExonRefs().get(INTRON);

            if(transExonRefs == null)
                continue;

            // goal is to set the correct exon rank for these intronic reads based on the orientation of the splice junction
            for(TransExonRef transExonRef : transExonRefs)
            {

            }



            // mTransExonRefs.get(seIndex).putAll(read.getTransExonRefs());
        }

        for(final TranscriptData transData : transDataList)
        {
            if(!positionWithin(position, transData.TransStart, transData.TransEnd))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = transData.exons().get(i + 1);

                if(exon.ExonEnd < position && position < nextExon.ExonStart)
                {
                    int exonRank = orientation == 1 ? exon.ExonRank : nextExon.ExonRank;

                    TransExonRef transExonRef = new TransExonRef(transData.GeneId, transData.TransId, transData.TransName, exonRank);
                    mTransExonRefs.get(seIndex).put(INTRON, Lists.newArrayList(transExonRef));
                    break;
                }
            }
        }

        if(!mTransExonRefs.get(seIndex).isEmpty() && mRegionMatchTypes[seIndex] == RegionMatchType.NONE)
            mRegionMatchTypes[seIndex] = INTRON;
    }

    public void populateDiscordantTransExonRefs(int geneCollectionId, int seIndex)
    {
        for (ReadRecord read : mReads)
        {
            if (!read.Chromosome.equals(mChromosomes[seIndex]) || read.getGeneCollecton() != geneCollectionId)
                continue;

            mTransExonRefs.get(seIndex).putAll(read.getTransExonRefs());
        }
    }

    */

}
