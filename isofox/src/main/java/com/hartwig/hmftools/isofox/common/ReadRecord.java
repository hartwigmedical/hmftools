package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.generateMappedCoords;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MULTI_MAP_QUALITY_THRESHOLD;
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

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.bam.ClippedSide;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadRecord
{
    public final String Id;
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;

    public final String ReadBases;
    public final int Length; // of bases
    public final Cigar Cigar;

    private int mFlags;
    private String mMateChromosome;
    private int mMatePosStart;

    private final int[] mGeneCollections;
    private final boolean[] mIsGenicRegion;
    private final List<int[]> mMappedCoords;
    private boolean mLowerInferredAdded;
    private boolean mUpperInferredAdded;
    private final int[] mSoftClipRegionsMatched;
    private int mFragmentInsertSize;
    private String mSupplementaryAlignment;
    private boolean mHasInterGeneSplit;
    private short mMapQuality;
    private byte[] mBaseQualities;

    private int[] mJunctionPositions; // chimeric junctions

    private final Map<RegionReadData,RegionMatchType> mMappedRegions; // regions related to this read and their match type
    private final Map<Integer,TransMatchType> mTranscriptClassification;
    private final Map<RegionMatchType,List<TransExonRef>> mTransExonRefs;
    private final Map<RegionMatchType,List<TransExonRef>> mUpperTransExonRefs; // TE refs for upper coords if a spanning read

    public static final int NO_GENE_ID = -1;

    public static ReadRecord from(final SAMRecord record)
    {
        final String readId = record.isSecondaryAlignment() ? String.format("%s_%s",
                record.getReadName(), record.getAttribute("HI")) : record.getReadName();

        ReadRecord read = new ReadRecord(
                readId, record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getCigar(), record.getInferredInsertSize(), record.getFlags(),
                record.getMateReferenceName(), record.getMateAlignmentStart());

        read.setSuppAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
        read.setMapQuality((short)record.getMappingQuality());
        return read;
    }

    public ReadRecord(
            final String id, final String chromosome, int posStart, int posEnd, final String readBases, @NotNull final Cigar cigar,
            int insertSize, int flags, final String mateChromosome, int matePosStart)
    {
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Length = ReadBases.length();
        Cigar = cigar;

        mFlags = flags;
        mMateChromosome = mateChromosome;
        mMatePosStart = matePosStart;

        mGeneCollections = new int[] { NO_GENE_ID, NO_GENE_ID };
        mIsGenicRegion = new boolean[] { false, false };

        List<int[]> mappedCoords = generateMappedCoords(Cigar, PosStart);
        mMappedCoords = Lists.newArrayListWithCapacity(mappedCoords.size());
        mMappedCoords.addAll(mappedCoords);

        mMappedRegions = Maps.newHashMap();
        mTransExonRefs = Maps.newHashMap();
        mUpperTransExonRefs = Maps.newHashMap();
        mTranscriptClassification = Maps.newHashMap();
        mLowerInferredAdded = false;
        mUpperInferredAdded = false;
        mSoftClipRegionsMatched = new int[] {0, 0};
        mFragmentInsertSize = insertSize;
        mSupplementaryAlignment = null;
        mHasInterGeneSplit = false;
        mMapQuality = 0;
        mJunctionPositions = null;
        mBaseQualities = null;
    }

    public int range() { return PosEnd - PosStart; }

    public byte orientation() { return !isReadReversed() ? POS_ORIENT : NEG_ORIENT; }

    public int flags() { return mFlags; }
    public boolean isReadPaired() { return (mFlags & SAMFlag.READ_PAIRED.intValue()) != 0; }
    public boolean isReadReversed() { return (mFlags & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mFlags & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isDuplicate() { return (mFlags & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean isTranslocation() { return !Chromosome.equals(mMateChromosome); }
    public boolean isMateNegStrand() { return (mFlags & SAMFlag.MATE_REVERSE_STRAND.intValue()) != 0; }
    public boolean isMateUnmapped() { return (mFlags & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }
    public boolean isInversion() { return isReadReversed() == isMateNegStrand(); }
    public boolean isProperPair() { return (mFlags & SAMFlag.PROPER_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mFlags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isSecondaryAlignment() { return (mFlags & SAMFlag.SECONDARY_ALIGNMENT.intValue()) != 0; }

    public void setFragmentInsertSize(int size) { mFragmentInsertSize = size; }
    public void setSuppAlignment(final String suppAlign) { mSupplementaryAlignment = suppAlign; }
    public String getSuppAlignment() { return mSupplementaryAlignment; }

    public String getSuppAlignmentCsv()
    {
        return mSupplementaryAlignment != null ? mSupplementaryAlignment.replaceAll(",", ITEM_DELIM) : "NONE";
    }

    public boolean hasSuppAlignment() { return mSupplementaryAlignment != null; }

    public static ClippedSide clippedSide(final ReadRecord read)
    {
        // considers hard or soft clips
        boolean leftClipped = read.getSoftClipRegionsMatched()[SE_START] == 0 && read.Cigar.isLeftClipped();
        boolean rightClipped = read.getSoftClipRegionsMatched()[SE_END] == 0 && read.Cigar.isRightClipped();
        return ClippedSide.from(read.Cigar, leftClipped, rightClipped);
    }

    public int[] getSoftClipRegionsMatched() { return mSoftClipRegionsMatched; }

    public boolean isSoftClipped(int se)
    {
        if(mSoftClipRegionsMatched[se] > 0)
            return false;

        return se == SE_START ? leftSoftClipped(Cigar) : rightSoftClipped(Cigar);
    }

    public boolean containsSoftClipping() { return Cigar.containsOperator(CigarOperator.S); }

    public void setMapQuality(short mapQuality) { mMapQuality = mapQuality; }
    public short mapQuality() { return mMapQuality; }

    public void setBaseQualities(final byte[] qualities) { mBaseQualities = qualities; }
    public byte[] baseQualities() { return mBaseQualities; }

    public boolean isMultiMapped() { return mMapQuality <= MULTI_MAP_QUALITY_THRESHOLD; }

    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public String mateChromosome() { return mMateChromosome; }
    public int mateStartPosition() { return mMatePosStart; }

    public int[] getGeneCollectons() { return mGeneCollections; }
    public final boolean[] getIsGenicRegion() { return mIsGenicRegion; }

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

    public boolean matches(final ReadRecord other)
    {
        return Id.equals(other.Id) && Cigar.toString().equals(other.Cigar.toString()) && PosStart == other.PosStart && PosEnd == other.PosEnd;
    }

    public boolean spansGeneCollections()
    {
        return mGeneCollections[SE_START] != mGeneCollections[SE_END];
    }

    public void setFlag(SAMFlag flag, boolean toggle)
    {
        if(toggle)
            mFlags |= flag.intValue();
        else
            mFlags &= ~flag.intValue();
    }

    public void setStrand(boolean readReversed, boolean mateReadReversed)
    {
        setFlag(SAMFlag.READ_REVERSE_STRAND, readReversed);
        setFlag(SAMFlag.MATE_REVERSE_STRAND, mateReadReversed);
    }

    public final Map<RegionMatchType,List<TransExonRef>> getReadTransExonRefs() { return mTransExonRefs; }
    public final Map<RegionMatchType,List<TransExonRef>> getReadTransExonRefs(int se)
    {
        if(spansGeneCollections())
            return se == SE_START ? mTransExonRefs : mUpperTransExonRefs;
        else
            return mTransExonRefs;
    }

    public boolean isChimeric()
    {
        if(isTranslocation() || isInversion())
            return true;

        if(!isProperPair() || isSupplementaryAlignment() || mSupplementaryAlignment != null)
            return true;

        return false;
    }

    public List<int[]> getMappedRegionCoords() { return mMappedCoords; }

    public List<int[]> getMappedRegionCoords(boolean includeInferred)
    {
        if(includeInferred || (!mLowerInferredAdded && !mUpperInferredAdded))
            return mMappedCoords;

        final List<int[]> regions = Lists.newArrayList(mMappedCoords);

        if(mLowerInferredAdded)
            regions.remove(0);

        if(mUpperInferredAdded)
            regions.remove(regions.size() - 1);

        return regions;
    }

    public boolean overlapsMappedReads(int posStart, int posEnd)
    {
        return mMappedCoords.stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x[SE_START], x[SE_END]));
    }

    public int getCoordsBoundary(int se)
    {
        return se == SE_START ? mMappedCoords.get(0)[SE_START] : mMappedCoords.get(mMappedCoords.size() - 1)[SE_END];
    }

    public void processOverlappingRegions(final List<RegionReadData> regions)
    {
        // process all regions for each transcript as a group to look for inconsistencies with the transcript definition
        Set<Integer> transcripts = Sets.newHashSet();

        for(RegionReadData region : regions)
        {
            for(final TransExonRef ref : region.getTransExonRefs())
            {
                transcripts.add(ref.TransId);
            }

            RegionMatchType matchType = setRegionMatchType(region);
            mMappedRegions.put(region, matchType);

            boolean checkMissedJunctions = matchType == EXON_INTRON || (Cigar.containsOperator(CigarOperator.S) && exonBoundary(matchType));

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

                if (matchType == RegionMatchType.NONE)
                {
                    // should never happen since implies this read didn't hit the region at all
                    transMatchType = ALT;
                }
                else if (matchType == RegionMatchType.EXON_INTRON)
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

                for (int regionIndex = 0; regionIndex < transRegions.size(); ++regionIndex)
                {
                    RegionReadData region = transRegions.get(regionIndex);

                    int exonRank = region.getExonRank(transId);
                    maxExonRank = max(maxExonRank, exonRank);
                    minExonRank = minExonRank == 0 ? exonRank : min(exonRank, minExonRank);

                    int mappingIndex = getRegionMappingIndex(region);
                    int adjustedMappingIndex = mappingIndex;

                    if(mLowerInferredAdded)
                        --adjustedMappingIndex;

                    if (adjustedMappingIndex < 0 || adjustedMappingIndex != regionIndex)
                    {
                        transMatchType = ALT;
                        break;
                    }

                    RegionMatchType matchType = mMappedRegions.get(region);

                    if (matchType == RegionMatchType.EXON_INTRON)
                    {
                        transMatchType = TransMatchType.UNSPLICED;
                        break;
                    }
                    else
                    {
                        final int[] readSection = mMappedCoords.get(mappingIndex);
                        int readStartPos = readSection[SE_START];
                        int readEndPos = readSection[SE_END];

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

                if (transMatchType == UNKNOWN)
                {
                    int expectedRegions = maxExonRank - minExonRank + 1;
                    if (transRegions.size() < expectedRegions)
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
            // is classified as alt
            if(validTranscriptType(transMatchType) && containsSoftClipping() && !likelyAdaperSoftClipping())
            {
                if(leftSoftClipped(Cigar) && mSoftClipRegionsMatched[SE_START] == 0)
                    transMatchType = ALT;
                else if(rightSoftClipped(Cigar) && mSoftClipRegionsMatched[SE_END] == 0)
                    transMatchType = ALT;
            }

            mTranscriptClassification.put(transId, transMatchType);
        }
    }

    public boolean likelyAdaperSoftClipping()
    {
        return mFragmentInsertSize < Length;
    }

    public static final List<RegionReadData> getUniqueValidRegion(final ReadRecord read1, final ReadRecord read2)
    {
        final List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        final List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                .filter(x -> validExonMatch(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        for(RegionReadData region : regions2)
        {
            if (!regions.contains(region))
                regions.add(region);
        }

        return regions;
    }

    public static boolean validTranscriptType(TransMatchType transType)
    {
        return transType == EXONIC || transType == SPLICE_JUNCTION;
    }

    public int getRegionMappingIndex(final RegionReadData region)
    {
        for(int i = 0; i < mMappedCoords.size(); ++i)
        {
            final int[] readSection = mMappedCoords.get(i);

            if(positionsOverlap(readSection[SE_START], readSection[SE_END], region.start(), region.end()))
                return i;
        }

        return -1;
    }

    private RegionMatchType setRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = getRegionMappingIndex(region);
        if (mappingIndex < 0)
            return RegionMatchType.NONE;

        RegionMatchType matchType = getRegionMatchType(region, mappingIndex);
        mMappedRegions.put(region, matchType);
        return matchType;
    }

    public RegionMatchType getRegionMatchType(final RegionReadData region)
    {
        int mappingIndex = getRegionMappingIndex(region);
        if (mappingIndex < 0)
            return RegionMatchType.NONE;

        return getRegionMatchType(region, mappingIndex);
    }

    private RegionMatchType getRegionMatchType(final RegionReadData region, int mappingIndex)
    {
        if(mappingIndex < 0 || mappingIndex >= mMappedCoords.size())
            return RegionMatchType.NONE;

        final int[] readSection = mMappedCoords.get(mappingIndex);
        int readStartPos = readSection[SE_START];
        int readEndPos = readSection[SE_END];

        if (readEndPos < region.start() || readStartPos > region.end())
            return RegionMatchType.NONE;

        if (readStartPos < region.start() || readEndPos > region.end())
            return RegionMatchType.EXON_INTRON;

        if (readStartPos > region.start() && readEndPos < region.end())
            return WITHIN_EXON;

        return EXON_BOUNDARY;
    }

    public static void markRegionBases(final List<int[]> readCoords, final RegionReadData region)
    {
        int[] regionBaseDepth = region.refBasesMatched();

        if(regionBaseDepth == null)
            return;

        for(final int[] readSection : readCoords)
        {
            int readStartPos = readSection[SE_START];
            int readEndPos = readSection[SE_END];

            if (readStartPos > region.end() || readEndPos < region.start())
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

    public boolean containsSplit()
    {
        return Cigar.containsOperator(CigarOperator.N);
    }

    public boolean hasInterGeneSplit() { return mHasInterGeneSplit; }
    public void setHasInterGeneSplit() { mHasInterGeneSplit = true; }

    private static final int MIN_SC_BASE_MATCH = 2;
    public static final int MAX_SC_BASE_MATCH = 10;

    private void checkMissedJunctions(final RegionReadData region)
    {
        if(mSupplementaryAlignment != null)
            return;

        // check for reads either soft-clipped or seemingly unspliced, where the extra bases can match with the next exon

        // check start of read
        int[] readSection = mLowerInferredAdded ? mMappedCoords.get(1) : mMappedCoords.get(0);
        int readStartPos = readSection[SE_START];
        int readEndPos = readSection[SE_END];

        // don't understand why any deleted bases were taken into account, no longer appears correct or relevant
        // int deletedLength = Cigar.getCigarElements().stream().filter(x -> x.getOperator() == D).mapToInt(x -> x.getLength()).sum();

        int extraBaseLength = 0;
        int scLength = 0;

        boolean hasRegionOverhang = region.start() > readStartPos && readEndPos > region.start()
                && region.start() - readStartPos <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = region.start() - readStartPos;
        }

        if(Cigar.getFirstCigarElement().getOperator() == CigarOperator.S && readStartPos <= region.start())
        {
            scLength = Cigar.getFirstCigarElement().getLength();
            extraBaseLength += scLength;
        }

        // less any deleted bases
        // extraBaseLength = max(extraBaseLength - deletedLength, 0);

        // allow a single base match if only 1 region matches
        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // first check for a match with the next exon on the lower side
            final String extraBases = ReadBases.substring(0, extraBaseLength);

            final List<RegionReadData> matchedRegions = region.getPreRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, false)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mSoftClipRegionsMatched[SE_START] = matchedRegions.size();
                mMappedRegions.put(region, EXON_BOUNDARY);

                if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    // truncate the read positions back to match the exon boundary
                    if (!mLowerInferredAdded && hasRegionOverhang)
                        readSection[SE_START] += region.start() - readStartPos;
                }

                // if only one region is matched or the min bases matched is satisfied, then create a mapping to the next region,
                // otherwise treat the splice support as ambiguous (it not mapped to the next region)
                if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for (RegionReadData preRegion : matchedRegions)
                    {
                        // add matched coordinates for this exon and add it as a region
                        mMappedRegions.put(preRegion, EXON_BOUNDARY);
                        addInferredMappingRegion(true, preRegion.end() - extraBaseLength + 1, preRegion.end());
                    }
                }
            }
        }

        // check end of read
        readSection = mUpperInferredAdded ? mMappedCoords.get(mMappedCoords.size() - 2) : mMappedCoords.get(mMappedCoords.size() - 1);
        readStartPos = readSection[SE_START];
        readEndPos = readSection[SE_END];

        extraBaseLength = 0;
        scLength = 0;

        hasRegionOverhang = readEndPos > region.end() && readStartPos < region.end() && readEndPos - region.end() <= MAX_SC_BASE_MATCH;

        if(hasRegionOverhang)
        {
            extraBaseLength = readEndPos - region.end();
        }

        if(Cigar.getLastCigarElement().getOperator() == CigarOperator.S && readEndPos >= region.end())
        {
            scLength = Cigar.getLastCigarElement().getLength();
            extraBaseLength += scLength;
        }

        // extraBaseLength = max(extraBaseLength - deletedLength, 0);

        if(extraBaseLength >= 1 && extraBaseLength <= MAX_SC_BASE_MATCH && scLength <= MAX_SC_BASE_MATCH)
        {
            // now check for a match to the next exon up
            final String extraBases = ReadBases.substring(Length - extraBaseLength, Length);

            final List<RegionReadData> matchedRegions = region.getPostRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, true)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mSoftClipRegionsMatched[SE_END] = matchedRegions.size();

                mMappedRegions.put(region, EXON_BOUNDARY);

                if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength < MIN_SC_BASE_MATCH))
                {
                    if (!mUpperInferredAdded && hasRegionOverhang)
                        readSection[SE_END] -= readEndPos - region.end();
                }

                if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_SC_BASE_MATCH))
                {
                    for (RegionReadData postRegion : matchedRegions)
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

        final String otherRegionBases = matchToStart ? otherRegion.refBases().substring(0, extraBases.length())
                : otherRegion.refBases().substring(otherRegionLength - extraBases.length(), otherRegionLength);

        return (otherRegionBases.equals(extraBases));
    }

    private int mappedRegionCount()
    {
        // discount any inferred regions
        return mMappedCoords.size() - (mLowerInferredAdded ? 1 : 0) - (mUpperInferredAdded ? 1 : 0);
    }

    public boolean inferredCoordAdded(boolean isLower) { return isLower ? mLowerInferredAdded : mUpperInferredAdded; }

    private void addInferredMappingRegion(boolean isLower, int posStart, int posEnd)
    {
        if(isLower)
        {
            if (!mLowerInferredAdded)
            {
                mLowerInferredAdded = true;
                mMappedCoords.add(0, new int[] { posStart, posEnd });
            }
            else
            {
                // lengthen the new region if required
                int[] newSection = mMappedCoords.get(0);
                newSection[SE_START] = min(newSection[SE_START], posStart);
            }
        }
        else
        {
            if(!mUpperInferredAdded)
            {
                mUpperInferredAdded = true;
                mMappedCoords.add(new int[] {posStart, posEnd});
            }
            else
            {
                int[] newSection = mMappedCoords.get(mMappedCoords.size() - 1);
                newSection[SE_END] = max(newSection[SE_END], posEnd);
            }
        }
    }

    public final Map<RegionReadData,RegionMatchType> getMappedRegions() { return mMappedRegions; }

    public static List<RegionReadData> findOverlappingRegions(final List<RegionReadData> regions, final ReadRecord read)
    {
        return regions.stream()
                .filter(x -> read.overlapsMappedReads(x.start(), x.end()))
                .collect(Collectors.toList());
    }

    public void addIntronicTranscriptRefs(final List<TranscriptData> transDataList)
    {
        final List<TransExonRef> transRefList = Lists.newArrayList();

        for(final TranscriptData transData : transDataList)
        {
            if(!mMappedCoords.stream().anyMatch(x -> positionsWithin(x[SE_START], x[SE_END], transData.TransStart, transData.TransEnd)))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = transData.exons().get(i + 1);

                if(mMappedCoords.stream().anyMatch(x -> positionsWithin(x[SE_START], x[SE_END], exon.End, nextExon.Start)))
                {
                    int minExonRank = min(exon.Rank, nextExon.Rank);
                    transRefList.add(new TransExonRef(transData.GeneId, transData.TransId, transData.TransName, minExonRank));
                    break;
                }
            }
        }

        if(!transRefList.isEmpty())
            mTransExonRefs.put(INTRON, transRefList);
    }

    public final List<TransExonRef> getJunctionMatchingTransRefs(int junctionPosition, boolean isJunctionStart)
    {
        final List<TransExonRef> matchedTransRefs = Lists.newArrayList();

        mMappedRegions.entrySet().stream()
                .filter(x -> exonBoundary(x.getValue()))
                .filter(x -> (isJunctionStart && x.getKey().end() == junctionPosition)
                        || (!isJunctionStart && x.getKey().start() == junctionPosition))
                .forEach(x -> matchedTransRefs.addAll(x.getKey().getTransExonRefs()));

        return matchedTransRefs;
    }

    public final Map<Integer,TransMatchType> getTranscriptClassifications() { return mTranscriptClassification; }

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

    public String toString()
    {
        return String.format("range(%s: %d -> %d, range=%d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, range(), Length, Cigar != null ? Cigar.toString() : "", Id);
    }
}
