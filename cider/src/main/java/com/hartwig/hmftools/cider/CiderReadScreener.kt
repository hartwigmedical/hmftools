package com.hartwig.hmftools.cider;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.eclipse.collections.api.collection.ImmutableCollection;
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
    private CiderGeneDatastore mCiderGeneDatastore;

    private IAnchorBlosumSearcher mAnchorBlosumSearcher;
    int mMaxAnchorAlignDistance;
    int mMinAnchorOverlap;

    private Multimap<ReadKey, VJReadCandidate> mCdr3ReadMatchMap = ArrayListMultimap.create();

    // want a map to make sure we do not process same record twice
    private final Set<VJReadRecordKey> mProcessedReadRecords = new HashSet<>();

    private final List<SAMRecord> mAllMatchedReads = new ArrayList<>();

    public Cdr3ReadScreener(CiderGeneDatastore ciderGeneDatastore, IAnchorBlosumSearcher anchorBlosumSearcher,
            int maxAnchorAlignDistance, int minAnchorOverlap)
    {
        mCiderGeneDatastore = ciderGeneDatastore;
        mAnchorBlosumSearcher = anchorBlosumSearcher;
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
        else
        {
            // next match by direct compare, we always have to do both, reason is that if this is a
            // CDR3 gene then a read could match V, D, or J
            matchFound |= tryMatchByBlosum(samRecord);
        }

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

        for (VJAnchorReferenceLocation anchorLocation : mCiderGeneDatastore.getVJAnchorReferenceLocations())
        {
            @Nullable VJReadCandidate readCandidate = matchesAnchorLocation(samRecord, mapped, matchedGenes, anchorLocation);

            if (readCandidate != null)
            {
                anchorFound = true;
            }
        }

        if (!anchorFound)
        {
            int leftSoftClip = CigarUtils.leftSoftClip(samRecord);
            int rightSoftClip = CigarUtils.rightSoftClip(samRecord);

            if (leftSoftClip != 0 || rightSoftClip != 0)
            {
                for (IgTcrConstantRegion igTcrConstantRegion : mCiderGeneDatastore.getIgConstantRegions())
                {
                    // now try to match around location of constant regions
                    @Nullable VJReadCandidate readCandidate = tryMatchFromConstantRegion(samRecord, mapped, igTcrConstantRegion);

                    if (readCandidate != null)
                    {
                        anchorFound = true;
                    }
                }
            }
        }

        return anchorFound;
    }

    @Nullable
    public VJReadCandidate matchesAnchorLocation(final SAMRecord samRecord, final GenomeRegion mapped, final Set<VJGeneType> matchedGenes,
            final VJAnchorReferenceLocation anchorLocation)
    {
        int readLength = samRecord.getReadLength();

        if (!isRelevantToAnchorLocation(readLength, mapped, anchorLocation))
            return null;

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

            // want to make sure same gene is not included twice
            ImmutableCollection<VJAnchorTemplate> genes = mCiderGeneDatastore.getByAnchorGeneLocation(anchorLocation)
                    .select(o -> !matchedGenes.contains(o.getType()))
                    .toImmutableList();

            if (!genes.isEmpty())
            {
                matchedGenes.addAll(genes.stream().map(VJAnchorTemplate::getType).collect(Collectors.toList()));

                return createCdr3ReadMatch(samRecord,
                        genes,
                        VJReadCandidate.AnchorMatchMethod.ALIGN,
                        anchorLocation.getStrand() == Strand.REVERSE,
                        readAnchorStart,
                        readAnchorEnd, anchorLocation.getGeneLocation());
            }
        }
        return null;
    }

    private boolean tryMatchByBlosum(final SAMRecord samRecord)
    {
        AnchorBlosumMatch anchorBlosumMatch = null;

        for (Strand strand : Strand.values())
        {
            String readString = samRecord.getReadString();

            if (strand == Strand.REVERSE)
                readString = SequenceUtil.reverseComplement(readString);

            for (VJGeneType vjGeneType : VJGeneType.values())
            {
                anchorBlosumMatch = mAnchorBlosumSearcher.searchForAnchor(readString,
                        vjGeneType,
                        0,
                        samRecord.getReadLength());

                if (anchorBlosumMatch != null && anchorBlosumMatch.getSimilarityScore() > 0)
                {
                    VJReadCandidate readCandidate = createCdr3ReadMatch(samRecord,
                            anchorBlosumMatch.getTemplateGenes(),
                            VJReadCandidate.AnchorMatchMethod.BLOSUM,
                            strand == Strand.REVERSE,
                            anchorBlosumMatch.getAnchorStart(),
                            anchorBlosumMatch.getAnchorEnd(), null);

                    if (readCandidate != null)
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    /*
    public boolean tryMatchByExact(SAMRecord samRecord, Set<VJGeneType> matchedGenes)
    {
        String seq = samRecord.getReadString();
        String seqRevComp = SequenceUtil.reverseComplement(seq);
        boolean useRevComp = false;
        boolean matchFound = false;

        // for each record we need to find the anchor sequence and see if it matches
        for (String anchorSeq : mCiderGeneDatastore.getAnchorSequenceSet())
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
                ImmutableCollection<VJAnchorTemplate> genes = mCiderGeneDatastore.getByAnchorSequence(anchorSeq)
                        .select(o -> !matchedGenes.contains(o.getType()))
                        .toImmutableList();

                if (!genes.isEmpty())
                {
                    createCdr3ReadMatch(samRecord,
                            genes,
                            VJReadCandidate.AnchorMatchMethod.EXACT,
                            useRevComp,
                            anchorIndex,
                            anchorIndex + anchorSeq.length(),
                            null);

                    matchFound = true;
                }
            }
        }

        return matchFound;
    }*/

    @Nullable
    public VJReadCandidate tryMatchFromConstantRegion(
            final SAMRecord samRecord, final GenomeRegion mapped, IgTcrConstantRegion igTcrConstantRegion)
    {
        int readLength = samRecord.getReadLength();
        GeneLocation igcLocation = igTcrConstantRegion.getGeneLocation();

        // see if the anchor location is mapped around here
        if ((igcLocation.getPosStart() - readLength < mapped.end() &&
                igcLocation.getPosEnd() + readLength > mapped.start()) &&
                igcLocation.getChromosome().equals(mapped.chromosome()))
        {
            // this is mapped to an IG constant region
            // we want to then see if there is any clipping going on
            // only with soft clip can we say this is actually potentially a J region
            // If these are in positive strand:
            // VVVVV-DDDD-JJ-CCCC
            // If the read is mapped to the constant CCCC region, then the left soft clip might
            // contain J
            // in negative strand, it would be right soft clip

            AnchorBlosumMatch anchorBlosumMatch = null;

            if (igcLocation.getStrand() == Strand.FORWARD)
            {
                int leftSoftClip = CigarUtils.leftSoftClip(samRecord);

                if (leftSoftClip == 0)
                    return null;

                // now try to find an anchor here
                anchorBlosumMatch = mAnchorBlosumSearcher.searchForAnchor(samRecord.getReadString(),
                        igTcrConstantRegion.getCorrespondingJ(),
                        0,
                        leftSoftClip);
            }
            else
            {
                int rightSoftClip = CigarUtils.rightSoftClip(samRecord);

                if (rightSoftClip == 0)
                    return null;

                String reverseCompSeq = SequenceUtil.reverseComplement(samRecord.getReadString());
                anchorBlosumMatch = mAnchorBlosumSearcher.searchForAnchor(reverseCompSeq,
                        igTcrConstantRegion.getCorrespondingJ(), 0, rightSoftClip);
            }

            if (anchorBlosumMatch != null && anchorBlosumMatch.getSimilarityScore() > 0)
            {
                return createCdr3ReadMatch(samRecord,
                        anchorBlosumMatch.getTemplateGenes(),
                        VJReadCandidate.AnchorMatchMethod.BLOSUM,
                        igcLocation.getStrand() == Strand.REVERSE,
                        anchorBlosumMatch.getAnchorStart(),
                        anchorBlosumMatch.getAnchorEnd(), null);
            }
        }
        return null;
    }

    @Nullable
    public VJReadCandidate createCdr3ReadMatch(SAMRecord samRecord, ImmutableCollection<VJAnchorTemplate> vjAnchorTemplates,
            VJReadCandidate.AnchorMatchMethod templateMatchMethod, boolean useRevComp,
            int readAnchorStart, int readAnchorEnd, @Nullable GeneLocation templateLocation)
    {
        if (vjAnchorTemplates.isEmpty())
            return null;

        VJAnchorTemplate vjAnchorTemplate = vjAnchorTemplates.stream().findFirst().get();

        // find out the imgt gene type. They should be the same type
        VJGeneType geneType = vjAnchorTemplate.getType();

        // check to make sure all the same
        if (vjAnchorTemplates.stream().anyMatch(o -> o.getType() != geneType))
        {
            sLogger.error("multiple gene types found in same match: {}", vjAnchorTemplates);
            throw new RuntimeException("multiple gene types found in same match");
        }

        String templateAnchorAA = vjAnchorTemplate.getAnchorAminoAcidSequence();

        // since we don't actually know whether the aligned part is the anchor sequence, we have to use
        // the soft clip that we think make sense
        int leftSoftClip = CigarUtils.leftSoftClip(samRecord);
        int rightSoftClip = CigarUtils.rightSoftClip(samRecord);

        VJReadCandidate readMatch = new VJReadCandidate(
                samRecord, vjAnchorTemplates, geneType,
                vjAnchorTemplate.getAnchorSequence(),
                templateMatchMethod, useRevComp,
                readAnchorStart, readAnchorEnd, templateLocation, leftSoftClip, rightSoftClip);

        // set the similarity score
        readMatch.setSimilarityScore(BlosumSimilarityCalc.calcSimilarityScore(geneType.getVj(),
                        readMatch.getTemplateAnchorSequence(), readMatch.getAnchorSequence()));

        List<String> geneNames = vjAnchorTemplates.stream().map(o -> o.getName()).distinct().collect(Collectors.toList());

        sLogger.info("genes: {} read({}) match method({}) anchor range([{},{})) template loc({}) "
                        + "anchor AA({}) template AA({}) similarity({})",
                geneNames, samRecord, templateMatchMethod,
                readMatch.getAnchorOffsetStart(), readMatch.getAnchorOffsetEnd(),
                templateLocation,
                readMatch.getAnchorAA(), templateAnchorAA, readMatch.getSimilarityScore());

        // add it to list
        mCdr3ReadMatchMap.put(new ReadKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag()), readMatch);

        return readMatch;
    }

    public static boolean isRelevantToAnchorLocation(int readLength, final GenomeRegion mapped,
            final VJAnchorReferenceLocation anchorLocation)
    {
        if (!anchorLocation.getChromosome().equals(mapped.chromosome()))
            return false;

        // for V we only allow reads that are mapped upstream
        // for J we only allow reads that are mapped downstream
        // --------V-----------J--------
        // *****                  ******
        // reads mapped around the star sections are ok
        // we translate it to the genome coord space.

        int halfAnchorLength = anchorLocation.baseLength() / 2;
        boolean isMappedAroundHere;

        if (anchorLocation.getVj() == VJ.V && anchorLocation.getStrand() == Strand.FORWARD ||
                (anchorLocation.getVj() == VJ.J && anchorLocation.getStrand() == Strand.REVERSE))
        {
            // here we want the anchor to be downstream
            // we want anchor mapped to higher coord than or equal read
            // we allow anchor to overshoot the mapped region by half the anchor length
            //        |__read__|
            //    |_______________________________|     allowed anchor range
            // half anchor         read length
            isMappedAroundHere = anchorLocation.getStart() > mapped.start() - halfAnchorLength &&
                    anchorLocation.getEnd() < mapped.end() + readLength;
        }
        else
        {
            // here we want the anchor to be upstream
            // we want anchor mapped to lower coord than or equal read
            // we allow anchor to overshoot the mapped region by half the anchor length
            //                    |__read__|
            // |_______________________________|     allowed anchor range
            //   read length             half anchor
            isMappedAroundHere = anchorLocation.getEnd() < mapped.end() + halfAnchorLength &&
                    anchorLocation.getStart() > mapped.start() - readLength;
        }

        return isMappedAroundHere;
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
