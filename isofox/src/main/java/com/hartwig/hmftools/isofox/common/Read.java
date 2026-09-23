package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSize;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNKNOWN;

import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class Read
{
    private SAMRecord mRecord;

    // make private
    private int mPosStart;
    private int mPosEnd;

    private String mReadBases; // cached if trimmed
    private final String mOriginalCigarStr;
    private String mCigarStr;
    private final List<CigarElement> mCigarElements;

    public int mUnclippedStart;
    public int mUnclippedEnd;
    public final boolean mHasSplit;

    private final int[] mGeneCollections;
    private final boolean[] mIsGenicRegion;

    private final MappedCoords mMappedCoords;

    private SupplementaryReadData mSupplementaryData;
    private boolean mHasInterGeneSplit;
    private List<AltAlignment> mAltLoci; // alternate genomic mapping loci from XA tag; null if uniquely mapped
    private boolean mConsensusRead;

    private int[] mJunctionPositions; // chimeric junctions

    private final Map<RegionReadData,RegionMatchType> mMappedRegions; // regions related to this read and their match type
    private final Map<Integer,TransMatchType> mTranscriptClassification;
    private final Map<RegionMatchType,List<TransExonRef>> mTransExonRefs;

    public static final int NO_GENE_ID = -1;

    public Read(final SAMRecord record)
    {
        mRecord = record;
        mPosStart = record.getAlignmentStart();
        mPosEnd = record.getAlignmentEnd();

        mCigarElements = Lists.newArrayList(record.getCigar().getCigarElements());
        mOriginalCigarStr = record.getCigarString();
        mCigarStr = null; // set if trimmed

        mSupplementaryData = SupplementaryReadData.extractAlignment(mRecord);

        mHasSplit = mCigarElements.stream().anyMatch(x -> x.getOperator() == N);

        setBoundaries();

        mGeneCollections = new int[] { NO_GENE_ID, NO_GENE_ID };
        mIsGenicRegion = new boolean[] { false, false };

        mMappedCoords = MappedCoords.build(mCigarElements, mPosStart);

        mMappedRegions = Maps.newHashMap();
        mTransExonRefs = Maps.newHashMap();
        mTranscriptClassification = Maps.newHashMap();
        mHasInterGeneSplit = false;
        mJunctionPositions = null;

        mConsensusRead = record.hasAttribute(CONSENSUS_READ_ATTRIBUTE);
        mAltLoci = parseAltLoci(record.getStringAttribute(XA_ATTRIBUTE));
    }

    private void setBoundaries()
    {
        mUnclippedStart = mPosStart;

        if(!mCigarElements.isEmpty()&& mCigarElements.get(0).getOperator() == S)
            mUnclippedStart -= mCigarElements.get(0).getLength();

        mUnclippedEnd = mPosEnd;

        if(mCigarElements.size() >= 2)
        {
            int lastIndex = mCigarElements.size() - 1;
            if(mCigarElements.get(lastIndex).getOperator() == S)
                mUnclippedEnd += mCigarElements.get(lastIndex).getLength();
        }
    }

    public String id() { return mRecord.getReadName(); }
    public String chromosome() { return mRecord.getReferenceName(); }
    public int alignmentStart() { return mPosStart; }
    public int alignmentEnd() { return mPosEnd; }
    public int flags() { return mRecord.getFlags(); }

    public String mateChromosome() { return mRecord.getMateReferenceName(); }
    public int mateAlignmentStart() { return mRecord.getMateAlignmentStart(); }
    public int fragmentInsertSize() { return inferredInsertSize(mRecord); }
    public int mapQuality() { return mRecord.getMappingQuality(); }

    public byte orientByte() { return !isReadReversed() ? ORIENT_FWD : ORIENT_REV; }
    public Orientation orientation() { return !isReadReversed() ? Orientation.FORWARD : Orientation.REVERSE; }

    public List<CigarElement> cigarElements() { return mCigarElements; }
    public String cigarStr() { return mCigarStr != null ? mCigarStr : mOriginalCigarStr; }

    public String readBases() { return mReadBases != null ? mReadBases : mRecord.getReadString(); }

    public boolean containsSplit() { return mHasSplit; }

    public boolean isReadPaired() { return mRecord.getReadPairedFlag(); }
    public boolean isReadReversed() { return mRecord.getReadNegativeStrandFlag(); }
    public boolean isFirstOfPair() { return firstInPair(mRecord); }
    public boolean isDuplicate() { return mRecord.getDuplicateReadFlag(); }
    public boolean isTranslocation() { return !chromosome().equals(mateChromosome()); }
    public boolean isMateNegStrand() { return mateNegativeStrand(mRecord); }
    public boolean isMateUnmapped() { return mRecord.getMateUnmappedFlag(); }
    public boolean isInversion() { return isReadReversed() == isMateNegStrand(); }
    public boolean isSupplementaryAlignment() { return mRecord.getSupplementaryAlignmentFlag(); }

    public SupplementaryReadData supplementaryData() { return mSupplementaryData; }
    public boolean hasSuppAlignment() { return mSupplementaryData != null; }

    public MappedCoords mappedCoords() { return mMappedCoords; }

    // soft-clip methods which use the raw read's soft-clips, and not any exon-boundary inferred soft-clips
    public int unclippedStart() { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }
    public boolean isLeftClipped() { return mUnclippedStart != mPosStart; }
    public boolean isRightClipped() { return mUnclippedEnd != mPosEnd; }
    public boolean containsSoftClipping() { return isLeftClipped() || isRightClipped(); }
    public int leftClipLength() { return max(mPosStart - mUnclippedStart, 0); }
    public int rightClipLength() { return max(mUnclippedEnd - mPosEnd, 0); }
    public int longestSoftClip() { return max(leftClipLength(), rightClipLength()); }

    public Boolean longestSoftClipIsLeft()
    {
        int left = leftClipLength();
        int right = rightClipLength();

        if(left == 0 && right == 0)
            return null;

        return left >= right ? true : false;
    }

    public boolean isSoftClippedNoRegionMatch(int se)
    {
        if(mMappedCoords.isSoftClipRegionMatched(se))
            return false;

        return se == SE_START ? isLeftClipped() : isRightClipped();
    }

    public List<AltAlignment> altLoci() { return mAltLoci; }
    public int numLoci() { return mAltLoci != null ? 1 + mAltLoci.size() : 1; }

    public boolean isMultiMapped() { return numLoci() > 1; }

    public boolean isConsensusRead() { return mConsensusRead; }

    public int baseLength() { return readBases().length(); }

    public SAMRecord bamRecord() { return mRecord; }

    public int[] getGeneCollectons() { return mGeneCollections; }
    public boolean[] getIsGenicRegion() { return mIsGenicRegion; }

    public void setGeneCollection(int seIndex, int gc, boolean isGeneic)
    {
        mGeneCollections[seIndex] = gc;
        mIsGenicRegion[seIndex] = isGeneic;
    }

    public boolean withinGeneCollection() { return mIsGenicRegion[SE_START] && mIsGenicRegion[SE_END]; }
    public boolean overlapsGeneCollection() { return mIsGenicRegion[SE_START] || mIsGenicRegion[SE_END]; }

    public boolean fullyNonGenic()
    {
        return !mIsGenicRegion[SE_START] && !mIsGenicRegion[SE_END]
                && mGeneCollections[SE_START] != NO_GENE_ID && mGeneCollections[SE_END] != NO_GENE_ID;
    }

    public boolean matches(final Read other)
    {
        return id().equals(other.id()) && flags() == other.flags();
    }

    public boolean spansGeneCollections()
    {
        return mGeneCollections[SE_START] != mGeneCollections[SE_END];
    }

    public Map<RegionMatchType,List<TransExonRef>> getReadTransExonRefs() { return mTransExonRefs; }

    public boolean isChimeric()
    {
        if(isTranslocation() || isInversion())
            return true;

        if(isSupplementaryAlignment() || mSupplementaryData != null)
            return true;

        return false;
    }

    public List<BaseRegion> getMappedRegionCoords() { return mMappedCoords.alignments(); }
    public List<BaseRegion> getMappedRegionCoordsWithoutInferred() { return mMappedCoords.alignmentsWithoutInferred(); }

    public boolean overlapsMappedCoords(int posStart, int posEnd) { return mMappedCoords.alignmentsOverlap(posStart, posEnd); }
    public int getCoordsBoundary(int se) { return mMappedCoords.getCoordsBoundary(se); }

    // an alternate mapping locus from the XA tag: its genomic span and whether that alignment is spliced
    public static class AltAlignment
    {
        public final ChrBaseRegion Region;
        public final boolean Spliced;

        public AltAlignment(final ChrBaseRegion region, final boolean spliced)
        {
            Region = region;
            Spliced = spliced;
        }
    }

    private static List<AltAlignment> parseAltLoci(final String xaTag)
    {
        if(xaTag == null || xaTag.isEmpty())
            return null;

        List<AltAlignment> altLoci = Lists.newArrayList();

        // XA entry format: chr,(+/-)pos,CIGAR,NM
        for(String entry : xaTag.split(";"))
        {
            if(entry.isEmpty())
                continue;

            String[] fields = entry.split(",");

            if(fields.length < 2)
                continue;

            try
            {
                int position = Math.abs(Integer.parseInt(fields[1]));

                // alt CIGAR gives the ref span and splice status; fall back to a single base when absent
                int refLength = 1;
                boolean spliced = false;
                String cigarStr = fields.length > 2 ? fields[2] : null;

                if(cigarStr != null && !cigarStr.isEmpty())
                {
                    refLength = max(TextCigarCodec.decode(cigarStr).getReferenceLength(), 1);
                    spliced = cigarStr.indexOf('N') >= 0;
                }

                ChrBaseRegion region = new ChrBaseRegion(fields[0], position, position + refLength - 1);
                altLoci.add(new AltAlignment(region, spliced));
            }
            catch(NumberFormatException e)
            {
            }
        }

        return altLoci.isEmpty() ? null : altLoci;
    }

    public boolean likelyAdaperSoftClipping()
    {
        return fragmentInsertSize() < baseLength();
    }

    public RegionMatchType getRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = mMappedCoords.findRegionIndex(region);
        if(mappingIndex == MappedCoords.INVALID_INDEX)
            return RegionMatchType.NONE;

        return ReadTranscriptUtils.getRegionMatchType(mMappedCoords, region, mappingIndex);
    }

    public boolean hasInterGeneSplit() { return mHasInterGeneSplit; }
    public void setHasInterGeneSplit() { mHasInterGeneSplit = true; }

    public int mappedRegionCount() { return mMappedCoords.originalAlignmentCount(); } // does not include inferred regions

    public Map<RegionReadData,RegionMatchType> getMappedRegions() { return mMappedRegions; }

    public static List<RegionReadData> findOverlappingRegions(final List<RegionReadData> regions, final Read read)
    {
        return regions.stream()
                .filter(x -> read.overlapsMappedCoords(x.start(), x.end()))
                .collect(Collectors.toList());
    }

    public void addIntronicTranscriptRefs(final List<TranscriptData> transDataList)
    {
        List<TransExonRef> transRefList = Lists.newArrayList();

        for(TranscriptData transData : transDataList)
        {
            if(!mMappedCoords.alignmentsWithin(transData.TransStart, transData.TransEnd))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                ExonData exon = transData.exons().get(i);
                ExonData nextExon = transData.exons().get(i + 1);

                if(mMappedCoords.alignmentsWithin(exon.End, nextExon.Start))
                {
                    int minExonRank = min(exon.Rank, nextExon.Rank);
                    transRefList.add(new TransExonRef(
                            transData.GeneId, transData.TransId, transData.TransName, minExonRank, transData.IsCanonical));
                    break;
                }
            }
        }

        if(!transRefList.isEmpty())
            mTransExonRefs.put(INTRON, transRefList);
    }

    public List<TransExonRef> getJunctionMatchingTransRefs(int junctionPosition, boolean isJunctionStart)
    {
        List<TransExonRef> matchedTransRefs = Lists.newArrayList();

        mMappedRegions.entrySet().stream()
                .filter(x -> exonBoundary(x.getValue()))
                .filter(x -> (isJunctionStart && x.getKey().end() == junctionPosition)
                        || (!isJunctionStart && x.getKey().start() == junctionPosition))
                .forEach(x -> matchedTransRefs.addAll(x.getKey().getTransExonRefs()));

        return matchedTransRefs;
    }

    public Map<Integer,TransMatchType> getTranscriptClassifications() { return mTranscriptClassification; }

    public TransMatchType getTranscriptClassification(int transId)
    {
        TransMatchType transType = mTranscriptClassification.get(transId);
        return transType != null ? transType : UNKNOWN;
    }

    public int[] junctionPositions() { return mJunctionPositions; }

    public void setJunctionPosition(int se, int junctionPosition)
    {
        if(mJunctionPositions == null)
            mJunctionPositions = new int[SE_PAIR];

        mJunctionPositions[se] = junctionPosition;
    }

    public void trimAdapterSoftClipBases(int trimLength)
    {
        int softClipLength = orientation().isForward() ? rightClipLength() : leftClipLength();
        trimLength = min(trimLength, softClipLength);

        if(trimLength == 0)
            return;

        if(orientation().isForward())
        {
            // trim from upper end
            mReadBases = mRecord.getReadString().substring(0, mRecord.getReadBases().length - trimLength);

            mUnclippedEnd -= trimLength;

            int lastIndex = mCigarElements.size() - 1;

            if(trimLength < softClipLength)
            {
                mCigarElements.set(lastIndex, new CigarElement(softClipLength - trimLength, S));
            }
            else
            {
                mCigarElements.remove(lastIndex);
            }
        }
        else
        {
            // trim from lower end
            mReadBases = mRecord.getReadString().substring(trimLength);

            mUnclippedStart += trimLength;

            if(trimLength < softClipLength)
            {
                mCigarElements.set(0, new CigarElement(softClipLength - trimLength, S));
            }
            else
            {
                mCigarElements.remove(0);
            }
        }

        mCigarStr = cigarElementsToStr(mCigarElements);
    }

    public String toString()
    {
        return String.format("%s range(%s: %d -> %d) flags(%d) cigar(%s)",
                id(), chromosome(), mPosStart, mPosEnd, flags(), cigarStr());
    }

    @VisibleForTesting
    public void setFlag(SAMFlag flag, boolean toggle)
    {
        if(mRecord == null)
            return;

        int newFlags = mRecord.getFlags();

        if(toggle)
            newFlags |= flag.intValue();
        else
            newFlags &= ~flag.intValue();

        mRecord.setFlags(newFlags);
    }

    @VisibleForTesting
    public void setStrand(boolean readReversed, boolean mateReadReversed)
    {
        setFlag(SAMFlag.READ_REVERSE_STRAND, readReversed);
        setFlag(SAMFlag.MATE_REVERSE_STRAND, mateReadReversed);
    }

    @VisibleForTesting
    public void setSuppAlignment(final String suppAlign)
    {
        mRecord.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppAlign);
        mSupplementaryData = SupplementaryReadData.extractAlignment(mRecord);
    }
}
