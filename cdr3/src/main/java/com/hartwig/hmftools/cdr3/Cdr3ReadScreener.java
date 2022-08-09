package com.hartwig.hmftools.cdr3;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.samtools.CigarUtils;
import com.hartwig.hmftools.common.utils.IntPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.SequenceUtil;

public class Cdr3ReadScreener
{
    private static final Logger sLogger = LogManager.getLogger(Cdr3ReadScreener.class);

    // collect the reads and sort by types
    private VJGeneStore mVJGeneStore;
    int mMaxAnchorAlignDistance;
    int mMinAnchorOverlap;

    private Multimap<ReadKey, VJReadCandidate> mCdr3ReadMatchMap = ArrayListMultimap.create();

    // want a map to make sure we do not process same record twice
    private final Set<VJReadRecordKey> mProcessedReadRecords = new HashSet<>();

    private final List<SAMRecord> mAllMatchedReads = new ArrayList<>();

    public Cdr3ReadScreener(VJGeneStore vjGeneStore, int maxAnchorAlignDistance, int minAnchorOverlap)
    {
        mVJGeneStore = vjGeneStore;
        mMaxAnchorAlignDistance = maxAnchorAlignDistance;
        mMinAnchorOverlap = minAnchorOverlap;
    }

    public Multimap<ReadKey, VJReadCandidate> getVJReadCandidates()
    {
        return mCdr3ReadMatchMap;
    }

    public List<SAMRecord> getAllMatchedReads()
    {
        return mAllMatchedReads;
    }

    public void processSamRecord(SAMRecord samRecord)
    {
        // see if we already processed this read
        // note that we cannot use unmapped flag, we must check the reference index against -1
        // this is cause in many cases read record will use
        GenomeRegion mapped = samRecord.getReferenceIndex() == -1 ? null : GenomeRegions.create(samRecord.getReferenceName(),
                samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());

        VJReadRecordKey readRecordKey = new VJReadRecordKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag(), mapped);

        if (!mProcessedReadRecords.add(readRecordKey))
        {
            sLogger.trace("read record {} already processed", readRecordKey);
            return;
        }

        Set<VJGeneType> matchedGeneTypes = new HashSet<>();

        boolean matchFound = false;

        // first we try to match by genome region
        if (mapped != null)
        {
            matchFound |= tryMatchByAlignment(samRecord, mapped, matchedGeneTypes);
        }

        // next match by direct compare, we always have to do both, reason is that if this is a
        // CDR3 gene then a read could match V, D, or J
        matchFound |= tryMatchByExact(samRecord, matchedGeneTypes);

        if (matchFound)
        {
            mAllMatchedReads.add(samRecord);
        }
        else
        {
            sLogger.trace("read record {} does not match any VJ", samRecord);
        }
    }

    private boolean tryMatchByAlignment(final SAMRecord samRecord, final GenomeRegion mapped, Set<VJGeneType> matchedGenes)
    {
        boolean anchorFound = false;

        for (VJAnchorReferenceLocation anchorLocation : mVJGeneStore.getVJAnchorReferenceLocations())
        {
            @Nullable VJReadCandidate readCandidate = matchesAnchorLocation(samRecord, mapped, matchedGenes, anchorLocation);

            if (readCandidate != null)
            {
                mCdr3ReadMatchMap.put(new ReadKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag()), readCandidate);
                anchorFound = true;
            }
        }
        return anchorFound;
    }

    @Nullable
    public VJReadCandidate matchesAnchorLocation(final SAMRecord samRecord, final GenomeRegion mapped, final Set<VJGeneType> matchedGenes,
            final VJAnchorReferenceLocation anchorLocation)
    {
        int readLength = samRecord.getReadLength();

        // see if the anchor location is mapped around here
        if ((anchorLocation.getStart() - readLength < mapped.end() &&
            anchorLocation.getEnd() + readLength > mapped.start()) &&
            anchorLocation.getChromosome().equals(mapped.chromosome()))
        {
            int anchorLength = anchorLocation.baseLength();

            if (anchorLength != 30)
            {
                throw new RuntimeException("unexpected anchor length");
            }

            IntPair readAnchorRange = extrapolateAnchorReadRange(samRecord, anchorLocation);

            if (readAnchorRange == null)
            {
                return null;
            }

            int readAnchorStart = readAnchorRange.left;
            int readAnchorEnd = readAnchorRange.right;

            // as long as the record overlaps with the anchor location, we include it as candidate
            if (readAnchorStart + mMinAnchorOverlap < samRecord.getReadLength() && readAnchorEnd > mMinAnchorOverlap)
            {
                if (anchorLocation.getStrand() == Strand.REVERSE)
                {
                    // say if read length is 2, and read anchor start is 1, read anchor end is 2
                    // reverse it we want start = 1, end = 2
                    int revStart = samRecord.getReadLength() - readAnchorEnd;
                    readAnchorEnd = samRecord.getReadLength() - readAnchorStart;
                    readAnchorStart = revStart;
                }

                // we want to make sure same gene is not included twice
                List<VJGene> genes = mVJGeneStore.getByAnchorGeneLocation(anchorLocation).stream()
                        .filter(o -> !matchedGenes.contains(o))
                        .collect(Collectors.toList());

                if (!genes.isEmpty())
                {
                    matchedGenes.addAll(genes.stream().map(VJGene::getType).collect(Collectors.toList()));

                    return createCdr3ReadMatch(samRecord,
                            genes,
                            VJReadCandidate.AnchorMatchMethod.ALIGN,
                            anchorLocation.getStrand() == Strand.REVERSE,
                            readAnchorStart,
                            readAnchorEnd, anchorLocation.getGeneLocation());
                }
            }
        }
        return null;
    }

    public boolean tryMatchByExact(SAMRecord samRecord, Set<VJGeneType> matchedGenes)
    {
        String seq = samRecord.getReadString();
        String seqRevComp = SequenceUtil.reverseComplement(seq);
        boolean useRevComp = false;
        boolean matchFound = false;

        // for each record we need to find the anchor sequence and see if it matches
        for (String anchorSeq : mVJGeneStore.getAnchorSequenceSet())
        {
            int anchorIndex = seq.indexOf(anchorSeq);

            if (anchorIndex != -1)
            {
                useRevComp = false;
            }
            else if ((anchorIndex = seqRevComp.indexOf(anchorSeq)) != -1)
            {
                useRevComp = true;
            }

            if (anchorIndex != -1)
            {
                // want to make sure same gene is not included twice
                List<VJGene> genes = mVJGeneStore.getByAnchorSequence(anchorSeq).stream()
                        .filter(o -> !matchedGenes.contains(o.getType()))
                        .collect(Collectors.toList());

                if (!genes.isEmpty())
                {
                    VJReadCandidate readCandidate = createCdr3ReadMatch(samRecord,
                            genes,
                            VJReadCandidate.AnchorMatchMethod.EXACT,
                            useRevComp,
                            anchorIndex,
                            anchorIndex + anchorSeq.length(),
                            null);

                    mCdr3ReadMatchMap.put(new ReadKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag()), readCandidate);

                    matchFound = true;
                }
            }
        }

        return matchFound;
    }

    @Nullable
    public VJReadCandidate createCdr3ReadMatch(SAMRecord samRecord, List<VJGene> VJGenes,
            VJReadCandidate.AnchorMatchMethod templateMatchType, boolean useRevComp,
            int readAnchorStart, int readAnchorEnd, @Nullable GeneLocation templateLocation)
    {
        if (VJGenes.isEmpty())
            return null;

        String seq = samRecord.getReadString();

        if (useRevComp)
            seq = SequenceUtil.reverseComplement(seq);

        String anchorSeq = seq.substring(Math.max(readAnchorStart, 0), Math.min(readAnchorEnd, seq.length()));

        // find out the imgt gene type. They should be the same type
        VJGeneType geneType = VJGenes.stream().findFirst().get().getType();

        // check to make sure all the same
        if (VJGenes.stream().anyMatch(o -> o.getType() != geneType))
        {
            sLogger.error("multiple gene types found in same match: {}", VJGenes);
            throw new RuntimeException("multiple gene types found in same match");
        }

        // since we don't actually know whether the aligned part is the anchor sequence, we have to use
        // the soft clip that we think make sense
        int leftSoftClip = CigarUtils.leftSoftClip(samRecord);
        int rightSoftClip = CigarUtils.rightSoftClip(samRecord);

        VJReadCandidate readMatch = new VJReadCandidate(samRecord, VJGenes, geneType, templateMatchType, useRevComp,
                readAnchorStart, readAnchorEnd, templateLocation, leftSoftClip, rightSoftClip);

        List<String> geneNames = VJGenes.stream().map(o -> o.getName()).distinct().collect(Collectors.toList());

        sLogger.info("genes: {} read({}) match type({}) anchor range([{},{})) template loc({}) "
                        + "anchor AA({}) anchor seq({})",
                geneNames, samRecord, templateMatchType,
                readMatch.getAnchorOffsetStart(), readMatch.getAnchorOffsetEnd(),
                templateLocation,
                readMatch.getAnchorAA(), readMatch.getAnchorSequence());

        return readMatch;
    }

    // 0 based read offset
    // this is similar to SAMRecord getReadPositionAtReferencePosition
    // Two major differences are:
    // 1. this one returns 0 base offset
    // 2. this function will extrapolate into other regions if reference pos is not with any
    //    alignment block. For example, for cigar 20M1000N30M, if we query for a position
    //    10 bases after the 20M block, we will get a read offset of 30 (20M + 10 bases)
    static int extrapolateReadOffsetAtRefPosition(SAMRecord record, int referencePos)
    {
        // since it is almost impossible to have 3+ alignment blocks speed should not be an issue
        int closesAlignmentBlockDistance = Integer.MAX_VALUE;
        int readOffset = 0;

        for (AlignmentBlock alignmentBlock : record.getAlignmentBlocks())
        {
            int blockReferenceEnd = CoordMath.getEnd(alignmentBlock.getReferenceStart(), alignmentBlock.getLength());
            int blockReferenceStart = alignmentBlock.getReferenceStart();
            int blockDistance;

            if (referencePos < blockReferenceStart)
            {
                blockDistance = blockReferenceStart - referencePos;
            }
            else if (referencePos > blockReferenceEnd)
            {
                blockDistance = referencePos - blockReferenceStart;
            }
            else
            {
                // the anchor reference position is within this block
                blockDistance = 0;
            }

            if (blockDistance < closesAlignmentBlockDistance)
            {
                closesAlignmentBlockDistance = blockDistance;
                readOffset = alignmentBlock.getReadStart() + referencePos - blockReferenceStart - 1; // make it 0 based
            }
        }

        // if not within read then return -1
        if (readOffset < 0 || readOffset >= record.getReadLength())
        {
            return -1;
        }

        return readOffset;
    }

    @Nullable
    static IntPair extrapolateAnchorReadRange(SAMRecord record, final VJAnchorReferenceLocation anchorLocation)
    {
        // we always use the reference position to find it
        int anchorEndReferencePos = anchorLocation.anchorBoundaryReferencePosition();
        int anchorEndReadOffset = extrapolateReadOffsetAtRefPosition(record, anchorEndReferencePos);

        if (anchorEndReadOffset == -1)
        {
            return null;
        }
        else if (anchorEndReferencePos == anchorLocation.getStart())
        {
            return new IntPair(anchorEndReadOffset, anchorEndReadOffset + anchorLocation.baseLength());
        }
        else if (anchorEndReferencePos == anchorLocation.getEnd())
        {
            // make end exclusive
            return new IntPair(anchorEndReadOffset - anchorLocation.baseLength() + 1, anchorEndReadOffset + 1);
        }
        else
        {
            throw new IllegalStateException("anchor end reference pos must be either anchor start or anchor end");
        }
    }
}
