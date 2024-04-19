package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.FACING;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.common.IndelCoords;

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
        if(alignments.isEmpty())
            return;

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> zeroQualAlignments = Lists.newArrayList();

        for(AlignData alignment : alignments)
        {
            if(alignment.MapQual > 0)
                validAlignments.add(alignment);
            else
                zeroQualAlignments.add(alignment);
        }

        if(validAlignments.isEmpty())
            return;

        validAlignments.forEach(x -> x.setFullSequenceData(mAssemblyAlignment.fullSequence(), mAssemblyAlignment.fullSequenceLength()));

        if(validAlignments.size() == 1)
        {
            AlignData singleAlignment = validAlignments.get(0);

            // check for cigar-based INDELs and SGLs with soft-clips
            boolean formsIndel = formIndelBreakends(singleAlignment);

            if(!formsIndel && singleAlignment.maxSoftClipLength() < ALIGNMENT_MIN_SOFT_CLIP)
            {
                formSingleBreakend(validAlignments.get(0), zeroQualAlignments);
            }
        }
        else
        {
            // otherwise handle multiple alignments
            formMultipleBreakends(validAlignments, zeroQualAlignments);
        }
    }

    private boolean formIndelBreakends(final AlignData alignment)
    {
        // parse the CIGAR to get the indel coords
        IndelCoords indelCoords = findIndelCoords(alignment.RefLocation.start(), alignment.cigarElements(), MIN_INDEL_SUPPORT_LENGTH);

        if(indelCoords == null || indelCoords.Length < MIN_INDEL_LENGTH)
            return false;

        int indelSeqStart = alignment.sequenceStart() + indelCoords.PosStart - alignment.RefLocation.start();
        int indelSeqEnd = indelSeqStart + (indelCoords.isInsert() ? indelCoords.Length : 1);

        String insertedBases = "";

        if(indelCoords.isInsert())
        {
            insertedBases = mAssemblyAlignment.fullSequence().substring(indelSeqStart, indelSeqEnd);
        }

        HomologyData homology = null;

        Breakend lowerBreakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, indelCoords.PosStart, POS_ORIENT, insertedBases, homology);

        mAssemblyAlignment.addBreakend(lowerBreakend);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), indelSeqStart, POS_ORIENT, 1, alignment);

        lowerBreakend.addSegment(segment);

        Breakend upperBreakend = new Breakend(
                mAssemblyAlignment, alignment.RefLocation.Chromosome, indelCoords.PosEnd, NEG_ORIENT, insertedBases, homology);

        mAssemblyAlignment.addBreakend(upperBreakend);

        BreakendSegment nextSegment = new BreakendSegment(
                mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), indelSeqEnd, NEG_ORIENT, 1, alignment);

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
        byte orientation;
        int softClipLength;
        String insertedBases;

        if(alignment.leftSoftClipLength() >= alignment.rightSoftClipLength())
        {
            breakendPosition = alignment.RefLocation.start();
            orientation = NEG_ORIENT;
            softClipLength = alignment.leftSoftClipLength();
            insertedBases = fullSequence.substring(0, softClipLength);
        }
        else
        {
            breakendPosition = alignment.RefLocation.end();
            orientation = POS_ORIENT;
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
        boolean sqlRefEnd = alignment.orientation() == POS_ORIENT ? checkStart : !checkStart;

        int softClipLength = sqlRefEnd ? alignment.leftSoftClipLength() : alignment.rightSoftClipLength();

        if(softClipLength < ALIGNMENT_MIN_SOFT_CLIP)
            return false;

        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        int sglPosition = sqlRefEnd ? alignment.RefLocation.start() : alignment.RefLocation.end();

        String insertSequence = checkStart ?
                fullSequence.substring(0, softClipLength) : fullSequence.substring(fullSequenceLength - softClipLength);

        byte sglOrientation = segmentOrientation(alignment, !checkStart);

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

            byte breakendOrientation = segmentOrientation(alignment, true);
            int breakendPosition = alignment.isForward() ? alignment.RefLocation.end() : alignment.RefLocation.start();

            HomologyData homology = null;
            String insertedBases = "";

            AlignData nextAlignment = alignments.get(i + 1);

            if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
            {
                homology = HomologyData.determineHomology(alignment, nextAlignment, mRefGenome);
            }
            else if(alignment.sequenceEnd() < nextAlignment.sequenceStart() - 1)
            {
                int insertSeqStart = alignment.sequenceEnd() + 1;
                int insertSeqEnd = nextAlignment.sequenceStart();
                insertedBases = fullSequence.substring(insertSeqStart, insertSeqEnd);
            }

            Breakend breakend = new Breakend(
                    mAssemblyAlignment, alignment.RefLocation.Chromosome, breakendPosition, breakendOrientation, insertedBases, homology);

            mAssemblyAlignment.addBreakend(breakend);

            BreakendSegment segment = new BreakendSegment(
                    mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), alignment.sequenceEnd(),
                    breakendOrientation, nextSegmentIndex++, alignment);

            breakend.addSegment(segment);

            byte nextOrientation = segmentOrientation(nextAlignment, false);
            int nextPosition = nextAlignment.isForward() ? nextAlignment.RefLocation.start() : nextAlignment.RefLocation.end();

            Breakend nextBreakend = new Breakend(
                    mAssemblyAlignment, nextAlignment.RefLocation.Chromosome, nextPosition, nextOrientation, insertedBases, homology);

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

    protected static byte segmentOrientation(final AlignData alignment, boolean linksEnd)
    {
        return (linksEnd == (alignment.orientation() == POS_ORIENT)) ? POS_ORIENT : NEG_ORIENT;
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
        if(lowerBreakend.Orientation != NEG_ORIENT || upperBreakend.Orientation != POS_ORIENT)
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
