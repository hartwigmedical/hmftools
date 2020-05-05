package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;

import static htsjdk.samtools.CigarOperator.D;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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

    private int mGeneCollectionId;
    private final List<int[]> mMappedCoords;
    private boolean mLowerInferredAdded;
    private boolean mUpperInferredAdded;
    private final int[] mSoftClipRegionsMatched;
    private int mFragmentInsertSize;
    private String mSupplementaryAlignment;

    private static final String SUPPLEMENTARY_ATTRIBUTE = "SA";

    private final Map<RegionReadData,RegionMatchType> mMappedRegions; // regions related to this read and their match type
    private final Map<Integer,TransMatchType> mTranscriptClassification;
    private final Map<RegionMatchType,List<TransExonRef>> mTransExonRefs;

    public static ReadRecord from(final SAMRecord record)
    {
        ReadRecord read = new ReadRecord(
                record.getReadName(), record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getCigar(), record.getInferredInsertSize(), record.getFlags(),
                record.getMateReferenceName(), record.getMateAlignmentStart());

        read.setSuppAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
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

        mGeneCollectionId = -1;

        List<int[]> mappedCoords = generateMappedCoords(Cigar, PosStart);
        mMappedCoords = Lists.newArrayListWithCapacity(mappedCoords.size());
        mMappedCoords.addAll(mappedCoords);

        mMappedRegions = Maps.newHashMap();
        mTransExonRefs = Maps.newHashMap();
        mTranscriptClassification = Maps.newHashMap();
        mLowerInferredAdded = false;
        mUpperInferredAdded = false;
        mSoftClipRegionsMatched = new int[] {0, 0};
        mFragmentInsertSize = insertSize;
        mSupplementaryAlignment = null;
    }

    public int range() { return PosEnd - PosStart; }

    public byte orientation()
    {
        // first in pair has orientation of +1 if not reversed, and vice versa for the second in the pair
        if(isFirstOfPair())
            return !isReadReversed() ? 1 : (byte)-1;
        else
            return isReadReversed() ? (byte)-1 : 1;
    }

    public int flags() { return mFlags; }
    public boolean isReadReversed() { return (mFlags & SAMFlag.READ_REVERSE_STRAND.intValue()) != 0; }
    public boolean isFirstOfPair() { return (mFlags & SAMFlag.FIRST_OF_PAIR.intValue()) != 0; }
    public boolean isDuplicate() { return (mFlags & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean isTranslocation() { return !Chromosome.equals(mMateChromosome); }
    public boolean isMateNegStrand() { return (mFlags & SAMFlag.MATE_REVERSE_STRAND.intValue()) != 0; }
    public boolean isMateUnmapped() { return (mFlags & SAMFlag.MATE_UNMAPPED.intValue()) != 0; }
    public boolean isInversion() { return isReadReversed() == isMateNegStrand(); }
    public boolean isProperPair() { return (mFlags & SAMFlag.PROPER_PAIR.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mFlags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isSecondaryAlignment() { return (mFlags & SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue()) != 0; }

    public void setFragmentInsertSize(int size) { mFragmentInsertSize = size; }
    public void setSuppAlignment(final String suppAlign) { mSupplementaryAlignment = suppAlign; }
    public String getSuppAlignment() { return mSupplementaryAlignment; }
    public int fragmentInsertSize() { return mFragmentInsertSize; }

    public String mateChromosome() { return mMateChromosome; }
    public int mateStartPosition() { return mMatePosStart; }

    public int getGeneCollecton() { return mGeneCollectionId; }

    public void setFlag(SAMFlag flag, boolean toggle)
    {
        if(toggle)
            mFlags |= flag.intValue();
        else
            mFlags &= ~flag.intValue();
    }

    public void setStrand(boolean isNeg, boolean mateIsNeg)
    {
        setFlag(SAMFlag.READ_REVERSE_STRAND, isNeg);
        setFlag(SAMFlag.MATE_REVERSE_STRAND, mateIsNeg);
    }

    public final Map<RegionMatchType,List<TransExonRef>> getTransExonRefs() { return mTransExonRefs; }

    public static boolean isTranslocation(@NotNull final SAMRecord record)
    {
        return !record.getReferenceName().equals(record.getMateReferenceName());
    }

    public static boolean isInversion(@NotNull final SAMRecord record)
    {
        return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
    }

    public boolean isChimeric()
    {
        if(isTranslocation() || isInversion())
            return true;

        if(!isProperPair() && isSecondaryAlignment() || isSupplementaryAlignment() || mSupplementaryAlignment != null)
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

    public static final List<int[]> generateMappedCoords(final Cigar cigar, int posStart)
    {
        final List<int[]> mappedCoords = Lists.newArrayList();

        // first establish whether the read is split across 2 distant regions, and if so which it maps to
        int posOffset = 0;
        boolean continueRegion = false;

        for(CigarElement element : cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.S)
            {
                // nothing to skip
            }
            else if(element.getOperator() == D)
            {
                posOffset += element.getLength();
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.I)
            {
                // nothing to skip
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                posOffset += element.getLength();
                continueRegion = false;
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                int readStartPos = posStart + posOffset;
                int readEndPos = readStartPos + element.getLength() - 1;

                if(continueRegion && !mappedCoords.isEmpty())
                {
                    int[] lastRegion = mappedCoords.get(mappedCoords.size() - 1);
                    lastRegion[SE_END] = readEndPos;
                }
                else
                {
                    mappedCoords.add(new int[] { readStartPos, readEndPos });
                }

                posOffset += element.getLength();
                continueRegion = false;
            }
        }

        return mappedCoords;
    }

    public void processOverlappingRegions(final List<RegionReadData> regions)
    {
        // process all regions for each transcript as a group to look for inconsistencies with the transcript definition
        List<Integer> transcripts = Lists.newArrayList();

        for(RegionReadData region : regions)
        {
            for(final TransExonRef ref : region.getTransExonRefs())
            {
                if (!transcripts.contains(ref.TransId))
                    transcripts.add(ref.TransId);
            }

            RegionMatchType matchType = getRegionMatchType(region);
            mMappedRegions.put(region, matchType);

            boolean checkMissedJunctions = matchType == EXON_INTRON || (Cigar.containsOperator(CigarOperator.S) && exonBoundary(matchType));

            if(checkMissedJunctions)
                checkMissedJunctions(region);
        }

        for(int transId : transcripts)
        {
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

                    if(mLowerInferredAdded)
                        --mappingIndex;

                    if (mappingIndex < 0 || mappingIndex != regionIndex)
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

            mTranscriptClassification.put(transId, transMatchType);
        }
    }

    public void captureGeneInfo(int geneCollectionId)
    {
        for(Map.Entry<RegionReadData,RegionMatchType> entry : mMappedRegions.entrySet())
        {
            List<TransExonRef> transRefList = mTransExonRefs.get(entry.getValue());

            if(transRefList == null)
            {
                mTransExonRefs.put(entry.getValue(), Lists.newArrayList(entry.getKey().getTransExonRefs()));
            }
            else
            {
                transRefList.addAll(entry.getKey().getTransExonRefs());
            }
        }

        mMappedRegions.clear();
        mTranscriptClassification.clear();
        mGeneCollectionId = geneCollectionId;
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

    public static boolean hasSkippedExons(final List<RegionReadData> regions, int transId, int longFragmentLimit)
    {
        int minExonRank = -1;
        int maxExonRank = 0;
        int regionCount = 0;

        int minRegionPos = -1;
        int maxRegionPos = 0;

        for(RegionReadData region : regions)
        {
            if(!region.hasTransId(transId))
                continue;

            ++regionCount;
            int exonRank = region.getExonRank(transId);

            maxExonRank = max(maxExonRank, exonRank);
            minExonRank = minExonRank == -1 ? exonRank : min(exonRank, exonRank);
            minRegionPos = minRegionPos == -1 ? region.start() : min(minRegionPos, region.start());
            maxRegionPos = max(maxRegionPos, region.end());
        }

        int expectedRegions = maxExonRank - minExonRank + 1;
        return regionCount < expectedRegions && maxRegionPos - minRegionPos > longFragmentLimit * 2;
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
            int readStartPos = readSection[SE_START];
            int readEndPos = readSection[SE_END];

            if (!(readStartPos > region.end() || readEndPos < region.start()))
                return i;
        }

        return -1;
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

        if (readStartPos == region.start() && readEndPos == region.end())
            return EXON_MATCH;

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

    public static int calcFragmentLength(final ReadRecord read1, final ReadRecord read2)
    {
        int insertSize = abs(read1.fragmentInsertSize());

        if(!read1.containsSplit() && !read2.containsSplit())
            return insertSize;

        // find unique split lengths and remove them

        List<Integer> splitLengths = read1.Cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.N).map(x -> x.getLength()).collect(Collectors.toList());

        for(final CigarElement element : read2.Cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.N && !splitLengths.contains(element.getLength()))
                splitLengths.add(element.getLength());
        }

        int totalSplitLength = splitLengths.stream().mapToInt(x -> x).sum();

        return insertSize - totalSplitLength;
    }

    public boolean containsSplit()
    {
        return Cigar.containsOperator(CigarOperator.N);
    }

    public boolean containsSoftClipping()
    {
        return Cigar.containsOperator(CigarOperator.S);
    }

    private static int MIN_BASE_MATCH = 2;

    private void checkMissedJunctions(final RegionReadData region)
    {
        // check for reads either soft-clipped or apparently unspliced, where the extra bases can match with the next exon

        // check start of read
        int[] readSection = mLowerInferredAdded ? mMappedCoords.get(1) : mMappedCoords.get(0);
        int readStartPos = readSection[SE_START];
        int readEndPos = readSection[SE_END];

        int deletedLength = Cigar.getCigarElements().stream().filter(x -> x.getOperator() == D).mapToInt(x -> x.getLength()).sum();

        int extraBaseLength = 0;

        if(region.start() > readStartPos && readEndPos > region.start())
        {
            extraBaseLength = (int)(region.start() - readStartPos);
        }

        if(Cigar.getFirstCigarElement().getOperator() == CigarOperator.S && readStartPos == region.start())
        {
            extraBaseLength += Cigar.getFirstCigarElement().getLength();
        }

        // less any deleted bases
        extraBaseLength = max(extraBaseLength - deletedLength, 0);

        // allow a single base match if only 1 region matches
        if(extraBaseLength >= 1 && extraBaseLength <= 10)
        {
            // first check for a match with the next exon on the lower side
            final String extraBases = ReadBases.substring(0, extraBaseLength);

            final List<RegionReadData> matchedRegions = region.getPreRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, false)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mSoftClipRegionsMatched[SE_START] = matchedRegions.size();
                mMappedRegions.put(region, EXON_BOUNDARY);

                if (matchedRegions.size() > 1 && extraBaseLength < MIN_BASE_MATCH && region.start() > readStartPos)
                {
                    // treat the splice support as ambiguous and truncate the read positions
                    if (!mLowerInferredAdded)
                    {
                        readSection[SE_START] += region.start() - readStartPos;
                    }
                }
                else if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_BASE_MATCH))
                {
                    for (RegionReadData preRegion : matchedRegions)
                    {
                        // add matched coordinates for this exon and it as a region
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

        if(readEndPos > region.end() && readStartPos < region.end())
        {
            extraBaseLength = (int)(readEndPos - region.end());
        }

        if(Cigar.getLastCigarElement().getOperator() == CigarOperator.S && readEndPos == region.end())
        {
            extraBaseLength += Cigar.getLastCigarElement().getLength();
        }

        extraBaseLength = max(extraBaseLength - deletedLength, 0);

        if(extraBaseLength >= 1 && extraBaseLength <= 10)
        {
            // now check for a match to the next exon up
            final String extraBases = ReadBases.substring(Length - extraBaseLength, Length);

            final List<RegionReadData> matchedRegions = region.getPostRegions().stream()
                    .filter(x -> matchesOtherRegionBases(extraBases, x, true)).collect(Collectors.toList());

            if(!matchedRegions.isEmpty())
            {
                mSoftClipRegionsMatched[SE_END] = matchedRegions.size();

                mMappedRegions.put(region, EXON_BOUNDARY);

                if (matchedRegions.size() > 1 && extraBaseLength < MIN_BASE_MATCH && readEndPos > region.end())
                {
                    if (!mUpperInferredAdded)
                    {
                        readSection[SE_END] -= readEndPos - region.end();
                    }
                }
                else if (matchedRegions.size() == 1 || (matchedRegions.size() > 1 && extraBaseLength >= MIN_BASE_MATCH))
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

    public boolean isSoftClipped(int se)
    {
        if(mSoftClipRegionsMatched[se] > 0)
            return false;

        return se == SE_START ? Cigar.isLeftClipped() : Cigar.isRightClipped();
    }

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
                .filter(x -> read.overlapsMappedReads(x.PosStart, x.PosEnd))
                .collect(Collectors.toList());
    }

    public void addIntronicTranscriptRefs(final List<TranscriptData> transDataList)
    {
        final List<TransExonRef> transRefList = Lists.newArrayList();

        for(final TranscriptData transData : transDataList)
        {
            if(!positionsWithin(PosStart, PosEnd, transData.TransStart, transData.TransEnd))
                continue;

            for(int i = 0; i < transData.exons().size() - 1; ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = transData.exons().get(i + 1);

                if(exon.ExonEnd < PosStart && PosEnd < nextExon.ExonStart)
                {
                    int minExonRank = min(exon.ExonRank, nextExon.ExonRank);
                    transRefList.add(new TransExonRef(transData.GeneId, transData.TransId, transData.TransName, minExonRank));
                    break;
                }
            }
        }

        if(!transRefList.isEmpty())
            mTransExonRefs.put(INTRON, transRefList);
    }


    public final Map<Integer,TransMatchType> getTranscriptClassifications() { return mTranscriptClassification; }

    public TransMatchType getTranscriptClassification(int transId)
    {
        TransMatchType transType = mTranscriptClassification.get(transId);
        return transType != null ? transType : UNKNOWN;
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d, range=%d) length(%d) cigar(%s) id(%s)",
                Chromosome, PosStart, PosEnd, range(), Length, Cigar != null ? Cigar.toString() : "", Id);
    }

}
