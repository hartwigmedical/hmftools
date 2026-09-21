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
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNKNOWN;

import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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

        parseAltLoci(record.getStringAttribute(XA_ATTRIBUTE));
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

    public void processOverlappingRegions(final List<RegionReadData> regions)
    {
        // process all regions for each transcript as a group to look for inconsistencies with the transcript definition
        Set<Integer> transcripts = Sets.newHashSet();
        boolean hasSoftClipping = isLeftClipped() || isRightClipped();

        for(RegionReadData region : regions)
        {
            for(TransExonRef ref : region.getTransExonRefs())
            {
                transcripts.add(ref.TransId);
            }

            RegionMatchType matchType = setRegionMatchType(region);
            mMappedRegions.put(region, matchType);

            boolean checkMissedJunctions = matchType == EXON_INTRON || (hasSoftClipping && exonBoundary(matchType));

            if(checkMissedJunctions)
                checkMissedJunctions(region);
        }

        for(int transId : transcripts)
        {
            // determine for each transcript whether the mapped regions support a spliced transcript, unspliced or alternate splicing
            TransMatchType transMatchType = UNKNOWN;

            List<RegionReadData> transRegions = regions.stream()
                    .filter(x -> x.getTransExonRefs().stream().anyMatch(y -> y.TransId == transId))
                    .collect(Collectors.toList());

            // if any reads cross and exon-intron boundary, then mark the transcript as unspliced

            if(transRegions.size() == 1 && mappedRegionCount() == 1)
            {
                // simple case of a single exon and read section
                RegionReadData region = transRegions.get(0);
                RegionMatchType matchType = mMappedRegions.get(region);

                if(matchType == RegionMatchType.NONE)
                {
                    // should never happen since implies this read didn't hit the region at all
                    transMatchType = ALT;
                }
                else if(matchType == RegionMatchType.EXON_INTRON)
                {
                    transMatchType = TransMatchType.UNSPLICED;
                }
            }
            else if(mappedRegionCount() > transRegions.size())
            {
                transMatchType = ALT;
            }
            else
            {
                int minExonRank = 0;
                int maxExonRank = 0;

                Collections.sort(transRegions);

                for(int regionIndex = 0; regionIndex < transRegions.size(); ++regionIndex)
                {
                    RegionReadData region = transRegions.get(regionIndex);

                    int exonRank = region.getExonRank(transId);
                    maxExonRank = max(maxExonRank, exonRank);
                    minExonRank = minExonRank == 0 ? exonRank : min(exonRank, minExonRank);

                    int mappingIndex = mMappedCoords.findRegionIndex(region);
                    int adjustedMappingIndex = mappingIndex;

                    if(mMappedCoords.inferredAlignmentAdded(SE_START))
                        --adjustedMappingIndex;

                    if(adjustedMappingIndex < 0 || adjustedMappingIndex != regionIndex)
                    {
                        transMatchType = ALT;
                        break;
                    }

                    RegionMatchType matchType = mMappedRegions.get(region);

                    if(matchType == RegionMatchType.EXON_INTRON)
                    {
                        transMatchType = TransMatchType.UNSPLICED;
                        break;
                    }
                    else
                    {
                        BaseRegion readSection = mMappedCoords.regionByIndex(mappingIndex);
                        int readStartPos = readSection.start();
                        int readEndPos = readSection.end();

                        boolean missStart = readStartPos > region.start();
                        boolean missEnd = readEndPos < region.end();
                        if(regionIndex == 0)
                        {
                            if(missEnd)
                            {
                                transMatchType = ALT;
                                break;
                            }
                        }
                        else if(regionIndex == transRegions.size() - 1)
                        {
                            if(missStart)
                            {
                                transMatchType = ALT;
                                break;
                            }
                        }
                        else if(missStart || missEnd)
                        {
                            transMatchType = ALT;
                            break;
                        }
                    }
                }

                if(transMatchType == UNKNOWN)
                {
                    int expectedRegions = maxExonRank - minExonRank + 1;
                    if(transRegions.size() < expectedRegions)
                        transMatchType = ALT;
                }
            }

            if(transMatchType == UNKNOWN)
            {
                if(transRegions.size() > 1)
                {
                    transMatchType = SPLICE_JUNCTION;
                }
                else
                {
                    transMatchType = EXONIC;
                }
            }

            // any read with soft-clipping which cannot be mapped to the next exon, other than for short likely adapter sequence reads,
            // is classified as alt, unless the clip is within an exon and a threshold
            if(validTranscriptType(transMatchType) && containsSoftClipping() && !likelyAdaperSoftClipping())
            {
                if(isLeftClipped() && !mMappedCoords.isSoftClipRegionMatched(SE_START) && !shortClipWithinExon(SE_START, transRegions))
                    transMatchType = ALT;
                else if(isRightClipped() && !mMappedCoords.isSoftClipRegionMatched(SE_END) && !shortClipWithinExon(SE_END, transRegions))
                    transMatchType = ALT;
            }

            mTranscriptClassification.put(transId, transMatchType);
        }
    }

    private boolean shortClipWithinExon(int se, final List<RegionReadData> transRegions)
    {
        int clipLength = se == SE_START ? leftClipLength() : rightClipLength();

        if(clipLength > MAX_SC_WITHIN_EXON_LENGTH)
            return false;

        return se == SE_START ?
                getCoordsBoundary(SE_START) - clipLength >= transRegions.stream().mapToInt(RegionReadData::start).min().orElse(0)
                : getCoordsBoundary(SE_END) + clipLength <= transRegions.stream().mapToInt(RegionReadData::end).max().orElse(0);
    }

    private boolean likelyAdaperSoftClipping()
    {
        return fragmentInsertSize() < baseLength();
    }

    public static final List<RegionReadData> getUniqueValidRegion(final Read read1, final Read read2)
    {
        List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        for(RegionReadData region : regions2)
        {
            if(!regions.contains(region))
                regions.add(region);
        }

        return regions;
    }

    public static boolean validTranscriptType(TransMatchType transType)
    {
        return transType == EXONIC || transType == SPLICE_JUNCTION;
    }

    private RegionMatchType setRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = mMappedCoords.findRegionIndex(region);
        if(mappingIndex == MappedCoords.INVALID_INDEX)
            return RegionMatchType.NONE;

        RegionMatchType matchType = getRegionMatchType(region, mappingIndex);
        mMappedRegions.put(region, matchType);
        return matchType;
    }

    public RegionMatchType getRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = mMappedCoords.findRegionIndex(region);
        if(mappingIndex == MappedCoords.INVALID_INDEX)
            return RegionMatchType.NONE;

        return getRegionMatchType(region, mappingIndex);
    }

    private RegionMatchType getRegionMatchType(final RegionReadData region, int mappingIndex)
    {
        if(mappingIndex == MappedCoords.INVALID_INDEX || mappingIndex >= mMappedCoords.alignmentCount())
            return RegionMatchType.NONE;

        BaseRegion readSection = mMappedCoords.regionByIndex(mappingIndex);
        int readStartPos = readSection.start();
        int readEndPos = readSection.end();

        if(readEndPos < region.start() || readStartPos > region.end())
            return RegionMatchType.NONE;

        if(readStartPos < region.start() || readEndPos > region.end())
            return RegionMatchType.EXON_INTRON;

        if(readStartPos > region.start() && readEndPos < region.end())
            return WITHIN_EXON;

        return EXON_BOUNDARY;
    }

    public static void markRegionBases(final List<BaseRegion> readCoords, final RegionReadData region)
    {
        int[] regionBaseDepth = region.refBasesMatched();

        if(regionBaseDepth == null)
            return;

        for(BaseRegion readSection : readCoords)
        {
            int readStartPos = readSection.start();
            int readEndPos = readSection.end();

            if(readStartPos > region.end() || readEndPos < region.start())
                continue;

            // process this overlap
            int regionBaseIndex = readStartPos > region.start() ? readStartPos - region.start() : 0;
            int overlap = min(readEndPos, region.end()) - max(readStartPos, region.start()) + 1;

            if(regionBaseIndex + overlap > regionBaseDepth.length)
            {
                ISF_LOGGER.error("region({}) read coords({} -> {}) regionBaseIndex({}) overlap({}) regionLength({})",
                        region, readStartPos, readEndPos, regionBaseIndex, overlap, regionBaseDepth.length);
                return;
            }

            for(int j = regionBaseIndex; j < regionBaseIndex + overlap; ++j)
            {
                ++regionBaseDepth[j];
            }
        }
    }

    public boolean hasInterGeneSplit() { return mHasInterGeneSplit; }
    public void setHasInterGeneSplit() { mHasInterGeneSplit = true; }

    private static final int MIN_SC_BASE_MATCH = 2;
    public static final int MAX_SC_BASE_MATCH = 10;
    private static final int MAX_SC_WITHIN_EXON_LENGTH = 2; // must stay below the realignment window (REALIGN_MIN_SOFT_CLIP_BASE_LENGTH)

    private void checkMissedJunctions(final RegionReadData region)
    {
        if(mSupplementaryData != null)
            return;

        // check for reads either soft-clipped or seemingly unspliced, where the extra bases can match with the next exon

        // check start of read
        BaseRegion readSection = mMappedCoords.lowestAlignment(false);
        int readStartPos = readSection.start();
        int readEndPos = readSection.end();

        int extraBaseLength = 0;
        int scLength = 0;

        boolean hasRegionOverhang = region.start() > readStartPos && readEndPos > region.start()
                && region.start() - readStartPos <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = region.start() - readStartPos;
        }

        if(isLeftClipped() && readStartPos <= region.start())
        {
            scLength = leftClipLength();
            extraBaseLength += scLength;
        }

        // less any deleted bases
        // extraBaseLength = max(extraBaseLength - deletedLength, 0);

        // allow a single base match if only 1 region matches
        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // first check for a match with the next exon on the lower side
            String extraBases = readBases().substring(0, extraBaseLength);

            List<RegionReadData> matchedRegions = region.getPreRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, false)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mMappedCoords.addSoftClipRegionMatched(true, matchedRegions.size());
                mMappedRegions.put(region, EXON_BOUNDARY);

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    // truncate the read positions back to match the exon boundary
                    if(!mMappedCoords.lowerInferredAlignmentAdded() && hasRegionOverhang)
                        readSection.setStart(readSection.start() + region.start() - readStartPos);
                }

                // if only one region is matched or the min bases matched is satisfied, then create a mapping to the next region,
                // otherwise treat the splice support as ambiguous (it not mapped to the next region)
                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for(RegionReadData preRegion : matchedRegions)
                    {
                        // add matched coordinates for this exon and add it as a region
                        mMappedRegions.put(preRegion, EXON_BOUNDARY);
                        addInferredMappingRegion(true, preRegion.end() - extraBaseLength + 1, preRegion.end());
                    }
                }
            }
        }

        // check end of read
        readSection = mMappedCoords.highestAlignment(false);
        readStartPos = readSection.start();
        readEndPos = readSection.end();

        extraBaseLength = 0;
        scLength = 0;

        hasRegionOverhang = readEndPos > region.end() && readStartPos < region.end() && readEndPos - region.end() <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = readEndPos - region.end();
        }

        if(isRightClipped() && readEndPos >= region.end())
        {
            scLength = rightClipLength();
            extraBaseLength += scLength;
        }

        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // now check for a match to the next exon up
            int readLength = baseLength();
            String extraBases = readBases().substring(readLength - extraBaseLength, readLength);

            List<RegionReadData> matchedRegions = region.getPostRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, true)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mMappedCoords.addSoftClipRegionMatched(false, matchedRegions.size());

                mMappedRegions.put(region, EXON_BOUNDARY);

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    if(!mMappedCoords.upperInferredAlignmentAdded() && hasRegionOverhang)
                        readSection.setEnd(readSection.end() - (readEndPos - region.end()));
                }

                if(matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for(RegionReadData postRegion : matchedRegions)
                    {
                        mMappedRegions.put(postRegion, EXON_BOUNDARY);
                        addInferredMappingRegion(false, postRegion.start(), postRegion.start() + extraBaseLength - 1);
                    }
                }
            }
        }
    }

    private static boolean matchesOtherRegionBases(final String extraBases, final RegionReadData otherRegion, boolean matchToStart)
    {
        int otherRegionLength = otherRegion.length();

        if(extraBases.length() > otherRegionLength)
            return false;

        String otherRegionBases = matchToStart ? otherRegion.refBases().substring(0, extraBases.length())
                : otherRegion.refBases().substring(otherRegionLength - extraBases.length(), otherRegionLength);

        return (otherRegionBases.equals(extraBases));
    }

    private int mappedRegionCount() { return mMappedCoords.originalAlignmentCount(); } // does not include inferred regions

    private void addInferredMappingRegion(boolean isLower, int posStart, int posEnd)
    {
        mMappedCoords.addInferredRegion(isLower, posStart, posEnd);
    }

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
