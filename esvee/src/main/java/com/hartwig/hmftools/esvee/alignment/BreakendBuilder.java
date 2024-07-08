package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.FACING;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.CigarElement;

public class BreakendBuilder
{
    private final RefGenomeInterface mRefGenome;
    private final AssemblyAlignment mAssemblyAlignment;

    public BreakendBuilder(final RefGenomeInterface refGenome, final AssemblyAlignment assemblyAlignment)
    {
        mRefGenome = refGenome;
        mAssemblyAlignment = assemblyAlignment;
    }

    public void formBreakends(final List<AlignData> alignments)
    {
        try
        {
            doFormBreakends(alignments);
        }
        catch(Exception e)
        {
            SV_LOGGER.error("assemblyAlign({}) error: {}", mAssemblyAlignment, e.toString());
            e.printStackTrace();

            for(AlignData alignData : alignments)
            {
                SV_LOGGER.debug("alignment({})", alignData);
            }
        }
    }

    private void doFormBreakends(final List<AlignData> alignments)
    {
        if(alignments.isEmpty())
            return;

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        filterAlignments(alignments, validAlignments, lowQualAlignments);

        if(validAlignments.isEmpty())
            return;

        if(validAlignments.size() == 1)
        {
            AlignData singleAlignment = validAlignments.get(0);

            // check for cigar-based INDELs and SGLs with soft-clips
            boolean formsIndel = formIndelBreakends(singleAlignment);

            if(!formsIndel && singleAlignment.maxSoftClipLength() >= ALIGNMENT_MIN_SOFT_CLIP)
            {
                formSingleBreakend(validAlignments.get(0), lowQualAlignments);
            }
        }
        else
        {
            // otherwise handle multiple alignments
            formMultipleBreakends(validAlignments, lowQualAlignments);
        }
    }

    @VisibleForTesting
    public void filterAlignments(
            final List<AlignData> alignments, final List<AlignData> validAlignments, final List<AlignData> lowQualAlignments)
    {
        // exclude alignments with zero qual
        final List<AlignData> nonZeroAlignments = Lists.newArrayList();

        for(AlignData alignment : alignments)
        {
            if(alignment.MapQual > 0)
                nonZeroAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }

        if(nonZeroAlignments.isEmpty())
            return;

        // for all the rest calculated an adjusted alignment score by subtracting overlap (inexact homology) and repeated bases from the score
        String fullSequence = mAssemblyAlignment.fullSequence();

        nonZeroAlignments.forEach(x -> x.setFullSequenceData(fullSequence, mAssemblyAlignment.fullSequenceLength()));

        // set modified map qual and then filtered low qual alignments
        Collections.sort(nonZeroAlignments, Comparator.comparingInt(x -> x.sequenceStart()));

        for(int i = 0; i < nonZeroAlignments.size(); ++i)
        {
            AlignData alignment = nonZeroAlignments.get(i);

            int overlapStart = 0;
            int overlapEnd = 0;

            if(i > 0)
            {
                AlignData prevAlignment = nonZeroAlignments.get(i - 1);
                if(prevAlignment.sequenceEnd() >= alignment.sequenceStart())
                {
                    overlapStart = prevAlignment.sequenceEnd() - alignment.sequenceStart() + 1;
                }
            }

            if(i < nonZeroAlignments.size() - 1)
            {
                AlignData nextAlignment = nonZeroAlignments.get(i + 1);

                if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
                {
                    overlapEnd = alignment.sequenceEnd() - nextAlignment.sequenceStart() + 1;
                }
            }

            alignment.setAdjustedAlignment(fullSequence, overlapStart, overlapEnd);
        }

        for(AlignData alignment : nonZeroAlignments)
        {
            if(alignment.calcModifiedMapQual() >= ALIGNMENT_MIN_MOD_MAP_QUAL)
                validAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }
    }

    private boolean formIndelBreakends(final AlignData alignment)
    {
        // parse the CIGAR to get the indel coords, first looking for a close match to what was expected
        IndelCoords indelCoords = null;

        StructuralVariantType svType = mAssemblyAlignment.svType();

        if(svType == DEL || svType == DUP || svType == INS)
        {
            int svLength = mAssemblyAlignment.svLength();
            CigarElement specificIndel = new CigarElement(svLength, svType == DEL ? D : I);
            indelCoords = findIndelCoords(alignment.RefLocation.start(), alignment.cigarElements(), specificIndel);
        }

        // if no match was found, just take the longest
        if(indelCoords == null)
            indelCoords = findIndelCoords(alignment.RefLocation.start(), alignment.cigarElements(), MIN_INDEL_SUPPORT_LENGTH);

        if(indelCoords == null || indelCoords.Length < MIN_INDEL_LENGTH)
            return false;

        int indelSeqStart = alignment.sequenceStart() + indelCoords.PosStart - alignment.RefLocation.start();
        int indelSeqEnd = indelSeqStart + (indelCoords.isInsert() ? indelCoords.Length : 1);

        String insertedBases = "";
        HomologyData homology = null;

        int indelPosStart = indelCoords.PosStart;
        int indelPosEnd = indelCoords.PosEnd;

        if(indelCoords.isInsert())
        {
            if(indelSeqStart >= 0 && indelSeqEnd < mAssemblyAlignment.fullSequenceLength())
                insertedBases = mAssemblyAlignment.fullSequence().substring(indelSeqStart, indelSeqEnd);
        }
        else
        {
            // check for exact homology at the bases either side of the delete
            int maxLength = indelCoords.Length - 1;
            int homPosStart = indelPosStart + 1;
            String basesStart = mRefGenome.getBaseString( alignment.RefLocation.Chromosome, homPosStart, homPosStart + maxLength);

            int homPosEnd = indelPosEnd;
            String basesEnd = mRefGenome.getBaseString(alignment.RefLocation.Chromosome, homPosEnd, homPosEnd + maxLength);

            homology = HomologyData.determineHomology(basesEnd, basesStart, basesEnd, maxLength);

            if(homology.Homology.isEmpty())
            {
                homology = null;
            }
            else
            {
                // in this case the delete does not include the overlapped homology bases, so both breakends need to be shifted foward
                // by the same amount, being the exact homology at the start
                // shift breakend positions forward by exact homology
                indelPosStart += abs(homology.ExactStart);
                indelPosEnd += abs(homology.ExactStart);
            }
        }

        Breakend lowerBreakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, indelPosStart, FORWARD, insertedBases, homology);

        mAssemblyAlignment.addBreakend(lowerBreakend);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), indelSeqStart, FORWARD, 1, alignment);

        lowerBreakend.addSegment(segment);

        Breakend upperBreakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, indelPosEnd, REVERSE, insertedBases, homology);

        mAssemblyAlignment.addBreakend(upperBreakend);

        BreakendSegment nextSegment = new BreakendSegment(
                mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), indelSeqEnd, REVERSE, 1, alignment);

        upperBreakend.addSegment(nextSegment);

        lowerBreakend.setOtherBreakend(upperBreakend);
        upperBreakend.setOtherBreakend(lowerBreakend);

        return true;
    }

    private void formSingleBreakend(final AlignData alignment, final List<AlignData> zeroQualAlignments)
    {
        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        int breakendPosition;
        Orientation orientation;
        int softClipLength;
        String insertedBases;

        if(alignment.leftSoftClipLength() >= alignment.rightSoftClipLength())
        {
            breakendPosition = alignment.RefLocation.start();
            orientation = REVERSE;
            softClipLength = alignment.leftSoftClipLength();
            insertedBases = fullSequence.substring(0, softClipLength);
        }
        else
        {
            breakendPosition = alignment.RefLocation.end();
            orientation = FORWARD;
            softClipLength = alignment.rightSoftClipLength();
            insertedBases = fullSequence.substring(fullSequenceLength - softClipLength);
        }

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, null);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), fullSequenceLength, alignment.sequenceStart(), orientation, 0, alignment);

        breakend.addSegment(segment);

        // check for alternative alignments
        if(!zeroQualAlignments.isEmpty())
        {
            List<AlternativeAlignment> altAlignments = Lists.newArrayList();
            zeroQualAlignments.forEach(x -> altAlignments.addAll(x.altAlignments()));
            breakend.setAlternativeAlignments(altAlignments);
        }

        mAssemblyAlignment.addBreakend(breakend);
    }

    private boolean checkOuterSingle(
            final AlignData alignment, boolean checkStart, int nextSegmentIndex, final List<AlignData> zeroQualAlignments)
    {
        // switch the reference coord and CIGAR end being check if the segment was reverse aligned
        boolean sqlRefEnd = alignment.orientation().isForward() ? checkStart : !checkStart;

        int softClipLength = sqlRefEnd ? alignment.leftSoftClipLength() : alignment.rightSoftClipLength();

        if(softClipLength < ALIGNMENT_MIN_SOFT_CLIP)
            return false;

        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        int sglPosition = sqlRefEnd ? alignment.RefLocation.start() : alignment.RefLocation.end();

        String insertSequence = checkStart ?
                fullSequence.substring(0, softClipLength) : fullSequence.substring(fullSequenceLength - softClipLength);

        Orientation sglOrientation = segmentOrientation(alignment, !checkStart);

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, sglPosition, sglOrientation, insertSequence, null);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), fullSequenceLength, alignment.sequenceStart(), sglOrientation, nextSegmentIndex, alignment);

        breakend.addSegment(segment);

        mAssemblyAlignment.addBreakend(breakend);

        List<AlternativeAlignment> altAlignments = checkStart ?
                buildAltAlignments(alignment, null, zeroQualAlignments) : buildAltAlignments(null, alignment, zeroQualAlignments);

        if(altAlignments != null && !altAlignments.isEmpty())
            breakend.setAlternativeAlignments(altAlignments);

        return true;
    }

    private void formMultipleBreakends(final List<AlignData> alignments, final List<AlignData> zeroQualAlignments)
    {
        // look for consecutive non-zero-MQ alignments from which to form breakends
        Collections.sort(alignments, Comparator.comparingInt(x -> x.sequenceStart()));

        String fullSequence = mAssemblyAlignment.fullSequence();

        int nextSegmentIndex = 0;

        // check for a single at the start or end of the alignments
        if(checkOuterSingle(alignments.get(0), true, nextSegmentIndex, zeroQualAlignments))
            nextSegmentIndex++;

        for(int i = 0; i < alignments.size() - 1; ++i)
        {
            AlignData alignment = alignments.get(i);

            Orientation breakendOrientation = segmentOrientation(alignment, true);
            int breakendPosition = alignment.isForward() ? alignment.RefLocation.end() : alignment.RefLocation.start();

            HomologyData homology = null;
            String insertedBases = "";

            AlignData nextAlignment = alignments.get(i + 1);

            if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
            {
                String assemblyOverlapBases = mAssemblyAlignment.overlapBases();

                if(assemblyOverlapBases.isEmpty())
                {
                    assemblyOverlapBases = fullSequence.substring(nextAlignment.sequenceStart(), alignment.sequenceEnd() + 1);
                }

                homology = HomologyData.determineHomology(assemblyOverlapBases, alignment, nextAlignment, mRefGenome);
            }
            else if(alignment.sequenceEnd() < nextAlignment.sequenceStart() - 1)
            {
                int insertSeqStart = alignment.sequenceEnd() + 1;
                int insertSeqEnd = nextAlignment.sequenceStart();
                insertedBases = fullSequence.substring(insertSeqStart, insertSeqEnd);
            }

            Orientation nextOrientation = segmentOrientation(nextAlignment, false);
            int nextPosition = nextAlignment.isForward() ? nextAlignment.RefLocation.start() : nextAlignment.RefLocation.end();

            HomologyData firstHomology = homology;
            HomologyData nextHomology = homology;

            if(homology != null)
            {
                boolean firstIsLower = alignment.isLowerAlignment(nextAlignment);
                boolean sameOrientation = breakendOrientation == nextOrientation;

                // correct to the convention where the lower breakend has the higher start values if they're not symmetrical
                if(!homology.isSymmetrical() && sameOrientation && !firstIsLower)
                {
                    firstHomology = homology.invert(true, false);
                }

                if(sameOrientation)
                {
                    boolean reversePositions = !homology.isSymmetrical() && firstIsLower;
                    boolean reverseBases = sameOrientation;
                    nextHomology = homology.invert(reversePositions, reverseBases);
                }

                // shift breakend positions where there is overlap, and move breakends back into the ref bases by the inexact homology
                breakendPosition += firstHomology.positionAdjustment(breakendOrientation);
                nextPosition += nextHomology.positionAdjustment(nextOrientation);
            }

            Breakend breakend = new Breakend(
                    mAssemblyAlignment, alignment.RefLocation.Chromosome, breakendPosition, breakendOrientation, insertedBases, firstHomology);

            mAssemblyAlignment.addBreakend(breakend);

            BreakendSegment segment = new BreakendSegment(
                    mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), alignment.sequenceEnd(),
                    breakendOrientation, nextSegmentIndex++, alignment);

            breakend.addSegment(segment);

            String nextInsertedBases = breakendOrientation != nextOrientation ?
                    insertedBases : Nucleotides.reverseComplementBases(insertedBases);

            Breakend nextBreakend = new Breakend(
                    mAssemblyAlignment, nextAlignment.RefLocation.Chromosome, nextPosition, nextOrientation, nextInsertedBases, nextHomology);

            mAssemblyAlignment.addBreakend(nextBreakend);

            BreakendSegment nextSegment = new BreakendSegment(
                    mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), nextAlignment.sequenceStart(),
                    nextOrientation, nextSegmentIndex++, nextAlignment);

            nextBreakend.addSegment(nextSegment);

            breakend.setOtherBreakend(nextBreakend);
            nextBreakend.setOtherBreakend(breakend);

            List<AlternativeAlignment> altAlignments = buildAltAlignments(alignment, nextAlignment, zeroQualAlignments);

            if(altAlignments != null && !altAlignments.isEmpty())
            {
                breakend.setAlternativeAlignments(altAlignments);
                nextBreakend.setAlternativeAlignments(altAlignments);
            }
        }

        checkOuterSingle(alignments.get(alignments.size() - 1), false, nextSegmentIndex, zeroQualAlignments);
    }

    protected static Orientation segmentOrientation(final AlignData alignment, boolean linksEnd)
    {
        return (linksEnd == (alignment.orientation().isForward())) ? FORWARD : REVERSE;
    }

    private static List<AlternativeAlignment> buildAltAlignments(
            final AlignData alignmentStart, AlignData alignmentEnd, final List<AlignData> zeroQualAlignments)
    {
        List<AlignData> relatedAlignments = null;

        if(alignmentStart != null && alignmentEnd != null)
        {
            relatedAlignments = zeroQualAlignments.stream()
                    .filter(x -> x.sequenceStart() > alignmentStart.sequenceStart() && x.sequenceStart() < alignmentEnd.sequenceStart())
                    .collect(Collectors.toList());
        }
        else if(alignmentStart != null)
        {
            relatedAlignments = zeroQualAlignments.stream()
                    .filter(x -> x.sequenceStart() < alignmentStart.sequenceStart())
                    .collect(Collectors.toList());
        }
        else
        {
            relatedAlignments = zeroQualAlignments.stream()
                    .filter(x -> x.sequenceStart() > alignmentEnd.sequenceStart())
                    .collect(Collectors.toList());
        }

        List<AlternativeAlignment> altAlignments = Lists.newArrayList();

        relatedAlignments.forEach(x -> altAlignments.addAll(x.altAlignments()));

        return altAlignments;
    }

    public static void formBreakendFacingLinks(final List<Breakend> breakends)
    {
        // may move this logic within breakend formation if assembly sequences cover chains of links
        for(int i = 0; i < breakends.size() - 1; ++i)
        {
            Breakend breakend = breakends.get(i);

            if(!breakend.passing())
                continue;

            for(int j = i + 1; j < breakends.size(); ++j)
            {
                Breakend nextBreakend = breakends.get(j);

                if(!nextBreakend.passing())
                    continue;

                if(nextBreakend.otherBreakend() == breakend) // to avoid DUP breakends also being considered facing
                    continue;

                if(!breakend.Chromosome.equals(nextBreakend.Chromosome) || nextBreakend.Position - breakend.Position > PHASED_ASSEMBLY_MAX_TI)
                    break;

                if(inFacingLink(breakend, nextBreakend))
                {
                    breakend.addFacingBreakend(nextBreakend);
                    nextBreakend.addFacingBreakend(breakend);
                }
            }
        }
    }

    private static boolean inFacingLink(final Breakend lowerBreakend, final Breakend upperBreakend)
    {
        if(lowerBreakend.Orient.isForward() || upperBreakend.Orient.isReverse())
            return false;

        PhaseSet lowerPhaseSet = lowerBreakend.assembly().assemblies().get(0).phaseSet();
        PhaseSet upperPhaseSet = upperBreakend.assembly().assemblies().get(0).phaseSet();

        if(lowerPhaseSet == null || lowerPhaseSet != upperPhaseSet)
            return false;

        // must have overlapping aligned segments
        boolean foundOverlap = false;

        for(BreakendSegment lowerSegment : lowerBreakend.segments())
        {
            for(BreakendSegment upperSegment : upperBreakend.segments())
            {
                if(lowerSegment.Alignment.RefLocation.overlaps(upperSegment.Alignment.RefLocation))
                {
                    foundOverlap = true;
                    break;
                }
            }

            if(foundOverlap)
                break;
        }

        if(!foundOverlap)
            return false;

        boolean foundMatchingLink = false;

        for(AssemblyLink assemblyLink : lowerPhaseSet.assemblyLinks())
        {
            if(assemblyLink.type() == FACING)
            {
                // facing link must overlap this region
                Junction firstJunction = assemblyLink.first().junction();
                Junction secondJunction = assemblyLink.second().junction();

                int lowerLinkPosition = min(firstJunction.Position, secondJunction.Position);
                int upperLinkPosition = max(firstJunction.Position, secondJunction.Position);

                if(positionsOverlap(lowerBreakend.Position, upperBreakend.Position, lowerLinkPosition, upperLinkPosition))
                {
                    foundMatchingLink = true;
                    break;
                }
            }
        }

        return foundMatchingLink;
    }
}
