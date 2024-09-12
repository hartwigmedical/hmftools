package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineSequenceCount;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.AssemblyUtils;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
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

        AlignmentFilters.filterAlignments(mAssemblyAlignment, alignments, validAlignments, lowQualAlignments);

        if(validAlignments.isEmpty())
            return;

        if(validAlignments.size() == 1)
        {
            AlignData singleAlignment = validAlignments.get(0);

            // check for cigar-based INDELs and SGLs with soft-clips
            boolean formsIndel = formIndelBreakends(singleAlignment);

            if(!formsIndel)
                formSingleBreakend(singleAlignment, lowQualAlignments);
        }
        else
        {
            // otherwise handle multiple alignments
            formMultipleBreakends(validAlignments, lowQualAlignments);
        }
    }

    private boolean formIndelBreakends(final AlignData alignment)
    {
        // parse the CIGAR to get the indel coords, first looking for a close match to what was expected
        IndelCoords indelCoords = null;

        boolean isChained = false;

        if(mAssemblyAlignment.phaseSet() != null)
        {
            if(mAssemblyAlignment.phaseSet().assemblyLinks().size() == 1)
            {
                AssemblyLink assemblyLink = mAssemblyAlignment.phaseSet().assemblyLinks().get(0);
                StructuralVariantType svType = assemblyLink.svType();

                if(svType == DEL || svType == DUP || svType == INS)
                {
                    int svLength = assemblyLink.length();
                    CigarElement specificIndel = new CigarElement(svLength, svType == DEL ? D : I);
                    indelCoords = IndelCoords.findMatchingIndelCoords(alignment.positionStart(), alignment.cigarElements(), specificIndel);
                }
            }
            else
            {
                isChained = true;
            }
        }

        // if no match was found, just take the longest
        if(indelCoords == null)
            indelCoords = findIndelCoords(alignment.positionStart(), alignment.cigarElements(), MIN_INDEL_LENGTH);

        if(indelCoords == null || indelCoords.Length < MIN_INDEL_LENGTH)
            return false;

        int indelSeqStart = alignment.sequenceStart() + indelCoords.PosStart - alignment.positionStart();
        int indelSeqEnd = indelSeqStart + (indelCoords.isInsert() ? indelCoords.Length : 1);

        if(!isChained)
        {
            // the indel must have sufficient bases either side of it to be called
            int leftAnchorLength = indelSeqStart + 1;
            int rightAnchorLength = alignment.segmentLength() - indelSeqEnd - 1;

            if(leftAnchorLength < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH || rightAnchorLength < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH)
                return false;
        }

        String insertedBases = "";
        HomologyData homology = null;

        int indelPosStart = indelCoords.PosStart;
        int indelPosEnd = indelCoords.PosEnd;

        if(indelCoords.isInsert())
        {
            if(indelSeqStart >= 0 && indelSeqEnd < mAssemblyAlignment.fullSequenceLength())
                insertedBases = mAssemblyAlignment.fullSequence().substring(indelSeqStart + 1, indelSeqEnd + 1);
        }
        else
        {
            // check for exact homology at the bases either side of the delete
            int maxLength = indelCoords.Length - 1;
            int homPosStart = indelPosStart + 1;
            String basesStart = mRefGenome.getBaseString( alignment.chromosome(), homPosStart, homPosStart + maxLength);

            int homPosEnd = indelPosEnd;
            String basesEnd = mRefGenome.getBaseString(alignment.chromosome(), homPosEnd, homPosEnd + maxLength);

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
                mAssemblyAlignment, alignment.chromosome(), indelPosStart, FORWARD, insertedBases, homology);

        mAssemblyAlignment.addBreakend(lowerBreakend);

        BreakendSegment segment = new BreakendSegment(mAssemblyAlignment.id(), indelSeqStart, FORWARD, 1, alignment);

        lowerBreakend.addSegment(segment);

        Breakend upperBreakend = new Breakend(
                mAssemblyAlignment, alignment.chromosome(), indelPosEnd, REVERSE, insertedBases, homology);

        mAssemblyAlignment.addBreakend(upperBreakend);

        BreakendSegment nextSegment = new BreakendSegment(mAssemblyAlignment.id(), indelSeqEnd, REVERSE, 1, alignment);

        upperBreakend.addSegment(nextSegment);

        lowerBreakend.setOtherBreakend(upperBreakend);
        upperBreakend.setOtherBreakend(lowerBreakend);

        return true;
    }

    private void formSingleBreakend(final AlignData alignment, final List<AlignData> zeroQualAlignments)
    {
        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        if(mAssemblyAlignment.assemblies().size() == 2 && mAssemblyAlignment.assemblies().stream().allMatch(x -> x.outcome() == LOCAL_INDEL))
        {
            // keep the breakends matching the original assembly and its local ref match
            processSglAsLocalIndel(alignment, zeroQualAlignments);
            return;
        }

        int breakendPosition;
        Orientation orientation;
        int softClipLength;
        String insertedBases;

        if(alignment.leftSoftClipLength() >= alignment.rightSoftClipLength())
        {
            breakendPosition = alignment.positionStart();
            orientation = REVERSE;
            softClipLength = alignment.leftSoftClipLength();
            insertedBases = fullSequence.substring(0, softClipLength);
        }
        else
        {
            breakendPosition = alignment.positionEnd();
            orientation = FORWARD;
            softClipLength = alignment.rightSoftClipLength();
            insertedBases = fullSequence.substring(fullSequenceLength - softClipLength);
        }

        int requiredSoftClipLength = isLineInsertion(insertedBases, orientation) ? LINE_MIN_EXTENSION_LENGTH : ALIGNMENT_MIN_SOFT_CLIP;

        if(softClipLength < requiredSoftClipLength)
            return;

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.chromosome(), breakendPosition, orientation, insertedBases, null);

        BreakendSegment segment = new BreakendSegment(mAssemblyAlignment.id(), alignment.sequenceStart(), orientation, 0, alignment);

        breakend.addSegment(segment);

        // check for alternative alignments
        if(!zeroQualAlignments.isEmpty())
        {
            List<AlternativeAlignment> altAlignments = Lists.newArrayList();
            zeroQualAlignments.forEach(x -> altAlignments.addAll(x.allAlignments()));
            breakend.setAlternativeAlignments(altAlignments);
        }

        mAssemblyAlignment.addBreakend(breakend);
    }

    private void processSglAsLocalIndel(final AlignData alignment, final List<AlignData> zeroQualAlignments)
    {
        // keep the breakends matching the original assembly and its local ref match
        JunctionAssembly posAssembly = mAssemblyAlignment.assemblies().stream().filter(x -> x.isForwardJunction()).findFirst().orElse(null);
        JunctionAssembly negAssembly = mAssemblyAlignment.assemblies().stream().filter(x -> x != posAssembly).findFirst().orElse(null);

        List<JunctionAssembly> assemblies = List.of(posAssembly, negAssembly);

        int sequenceStart = 0;
        int sequenceEnd = -1;
        int sequenceIndex = 0;
        for(JunctionAssembly assembly : assemblies)
        {
            boolean isOriginalAssembly = assembly.supportCount() > 0;

            Breakend breakend = new Breakend(
                    mAssemblyAlignment, assembly.junction().Chromosome, assembly.junction().Position,
                    assembly.junction().Orient, "", null);

            if(sequenceEnd > 0)
                sequenceStart = sequenceEnd + 1;

            sequenceEnd = sequenceStart + assembly.refBaseLength() - 1;

            int regionStart, regionEnd;

            if(assembly.isForwardJunction())
            {
                regionEnd = assembly.junction().Position;
                regionStart = isOriginalAssembly ? assembly.refBasePosition() : regionEnd - assembly.refBaseLength() + 1;
            }
            else
            {
                regionStart = assembly.junction().Position;
                regionEnd = isOriginalAssembly ? assembly.refBasePosition() : regionStart + assembly.refBaseLength() - 1;
            }

            ChrBaseRegion localRegion = new ChrBaseRegion(assembly.junction().Chromosome, regionStart, regionEnd);

            int score = isOriginalAssembly ? alignment.score() : regionEnd - regionStart + 1;

            AlignData breakendAlignment = new AlignData(
                    localRegion, sequenceStart, sequenceEnd, alignment.mapQual(), score, alignment.flags(), alignment.cigar(),
                    alignment.nMatches(), alignment.xaTag(), alignment.mdTag());

            breakendAlignment.setFullSequenceData(mAssemblyAlignment.fullSequence(), mAssemblyAlignment.fullSequenceLength());
            breakendAlignment.setAdjustedAlignment(mAssemblyAlignment.fullSequence(), 0, 0);

            BreakendSegment segment = new BreakendSegment(
                    mAssemblyAlignment.id(), sequenceStart, assembly.junction().Orient, sequenceIndex++, breakendAlignment);

            breakend.addSegment(segment);

            if(isOriginalAssembly && !zeroQualAlignments.isEmpty())
            {
                List<AlternativeAlignment> altAlignments = Lists.newArrayList();
                zeroQualAlignments.forEach(x -> altAlignments.addAll(x.allAlignments()));
                breakend.setAlternativeAlignments(altAlignments);
            }

            mAssemblyAlignment.addBreakend(breakend);
        }

        mAssemblyAlignment.breakends().get(0).setOtherBreakend(mAssemblyAlignment.breakends().get(1));
        mAssemblyAlignment.breakends().get(1).setOtherBreakend(mAssemblyAlignment.breakends().get(0));
    }

    private static boolean isLineInsertion(final String insertedBases, final Orientation orientation)
    {
        return findLineSequenceCount(
                insertedBases.getBytes(), 0, insertedBases.length() - 1,
                orientation.isForward() ? LINE_BASE_T : LINE_BASE_A) > 0;
    }

    private boolean checkOuterSingle(
            final AlignData alignment, boolean checkStart, int nextSegmentIndex, final List<AlignData> zeroQualAlignments)
    {
        // switch the reference coord and CIGAR end being check if the segment was reverse aligned
        boolean sglRefBaseAtEnd = alignment.orientation().isForward() ? checkStart : !checkStart;

        int softClipLength = sglRefBaseAtEnd ? alignment.leftSoftClipLength() : alignment.rightSoftClipLength();

        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        int sglPosition = sglRefBaseAtEnd ? alignment.positionStart() : alignment.positionEnd();

        String insertSequence = checkStart ?
                fullSequence.substring(0, softClipLength) : fullSequence.substring(fullSequenceLength - softClipLength);

        Orientation sglOrientation = segmentOrientation(alignment, !checkStart);

        int requiredSoftClipLength = isLineInsertion(insertSequence, sglOrientation) ? LINE_MIN_EXTENSION_LENGTH : ALIGNMENT_MIN_SOFT_CLIP;

        if(softClipLength < requiredSoftClipLength)
            return false;

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.chromosome(), sglPosition, sglOrientation, insertSequence, null);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), alignment.sequenceStart(), sglOrientation, nextSegmentIndex, alignment);

        breakend.addSegment(segment);

        mAssemblyAlignment.addBreakend(breakend);

        List<AlternativeAlignment> altAlignments = checkStart ?
                buildAltAlignments(alignment, null, zeroQualAlignments) : buildAltAlignments(null, alignment, zeroQualAlignments);

        if(altAlignments != null && !altAlignments.isEmpty())
            breakend.setAlternativeAlignments(altAlignments);

        return true;
    }

    private void checkAlignmentIndel(final AlignData alignment)
    {
        if(alignment.cigarElements().stream().anyMatch(x -> x.getOperator().isIndel() && x.getLength() >= MIN_INDEL_LENGTH))
            formIndelBreakends(alignment);
    }

    private void formMultipleBreakends(final List<AlignData> alignments, final List<AlignData> zeroQualAlignments)
    {
        // look for consecutive non-zero-MQ alignments from which to form breakends
        String fullSequence = mAssemblyAlignment.fullSequence();

        int nextSegmentIndex = 0;

        // check for a single at the start or end of the alignments
        if(checkOuterSingle(alignments.get(0), true, nextSegmentIndex, zeroQualAlignments))
            nextSegmentIndex++;

        for(int i = 0; i < alignments.size() - 1; ++i)
        {
            AlignData alignment = alignments.get(i);

            if(i == 0)
                checkAlignmentIndel(alignment);

            String breakendChr, nextChr;
            Orientation breakendOrientation, nextOrientation;
            int breakendPosition, nextPosition;

            if(alignment.hasSelectedAltAlignment())
            {
                AlternativeAlignment breakendAltAlignment = alignment.selectedAltAlignment();
                breakendChr = breakendAltAlignment.Chromosome;
                breakendPosition = breakendAltAlignment.Position;
                breakendOrientation = breakendAltAlignment.Orient;
            }
            else
            {
                breakendChr = alignment.chromosome();
                breakendPosition = alignment.isForward() ? alignment.positionEnd() : alignment.positionStart();
                breakendOrientation = segmentOrientation(alignment, true);
            }

            AlignData nextAlignment = alignments.get(i + 1);

            if(nextAlignment.hasSelectedAltAlignment())
            {
                AlternativeAlignment nextAltAlignment = nextAlignment.selectedAltAlignment();
                nextChr = nextAltAlignment.Chromosome;
                nextPosition = nextAltAlignment.Position;
                nextOrientation = nextAltAlignment.Orient;
            }
            else
            {
                nextChr = nextAlignment.chromosome();
                nextPosition = nextAlignment.isForward() ? nextAlignment.positionStart() : nextAlignment.positionEnd();
                nextOrientation = segmentOrientation(nextAlignment, false);
            }

            HomologyData homology = null;
            String insertedBases = "";

            if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
            {
                String assemblyOverlapBases = mAssemblyAlignment.overlapBases(alignment.sequenceEnd());

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

            String assemblyInsertedBases = breakendOrientation.isForward() ? insertedBases : Nucleotides.reverseComplementBases(insertedBases);

            Breakend breakend = new Breakend(
                    mAssemblyAlignment, breakendChr, breakendPosition, breakendOrientation, assemblyInsertedBases, firstHomology);

            mAssemblyAlignment.addBreakend(breakend);

            BreakendSegment segment = new BreakendSegment(
                    mAssemblyAlignment.id(), alignment.sequenceEnd(), breakendOrientation, nextSegmentIndex++, alignment);

            breakend.addSegment(segment);

            String nextInsertedBases = nextOrientation.isReverse() ? insertedBases : Nucleotides.reverseComplementBases(insertedBases);

            Breakend nextBreakend = new Breakend(
                    mAssemblyAlignment, nextChr, nextPosition, nextOrientation, nextInsertedBases, nextHomology);

            mAssemblyAlignment.addBreakend(nextBreakend);

            BreakendSegment nextSegment = new BreakendSegment(
                    mAssemblyAlignment.id(), nextAlignment.sequenceStart(), nextOrientation, nextSegmentIndex++, nextAlignment);

            nextBreakend.addSegment(nextSegment);

            breakend.setOtherBreakend(nextBreakend);
            nextBreakend.setOtherBreakend(breakend);

            List<AlternativeAlignment> altAlignments = buildAltAlignments(alignment, nextAlignment, zeroQualAlignments);

            if(altAlignments != null && !altAlignments.isEmpty())
            {
                breakend.setAlternativeAlignments(altAlignments);
                nextBreakend.setAlternativeAlignments(altAlignments);
            }

            checkAlignmentIndel(nextAlignment);
        }

        checkOuterSingle(alignments.get(alignments.size() - 1), false, nextSegmentIndex, zeroQualAlignments);

        // finally look for any facing breakends
        setFacingBreakends();
    }

    private void setFacingBreakends()
    {
        for(int i = 0; i < mAssemblyAlignment.breakends().size() - 1; ++i)
        {
            Breakend breakend = mAssemblyAlignment.breakends().get(i);
            Breakend nextBreakend = mAssemblyAlignment.breakends().get(i + 1);

            if(areFacingBreakends(breakend, nextBreakend))
            {
                breakend.addFacingBreakend(nextBreakend);
                nextBreakend.addFacingBreakend(breakend);
            }
        }

        // add links for any assemblies with multiple facing breakends in the phase set
        if(mAssemblyAlignment.phaseSet() == null || !mAssemblyAlignment.phaseSet().hasFacingLinks())
            return;

        for(AssemblyLink link : mAssemblyAlignment.phaseSet().assemblyLinks())
        {
            if(link.type() != LinkType.FACING)
                continue;

            Breakend firstBreakend = findAssemblyBreakend(link.first());
            Breakend secondBreakend = findAssemblyBreakend(link.second());

            if(firstBreakend != null && secondBreakend != null)
            {
                if(!firstBreakend.facingBreakends().contains(secondBreakend))
                    firstBreakend.addFacingBreakend(secondBreakend);

                if(!secondBreakend.facingBreakends().contains(firstBreakend))
                    secondBreakend.addFacingBreakend(firstBreakend);
            }
        }
    }

    private Breakend findAssemblyBreakend(final JunctionAssembly assembly)
    {
        Junction junction = assembly.junction();

        return mAssemblyAlignment.breakends().stream()
                .filter(x -> x.matchesCoordinates(junction.Chromosome, junction.Position, junction.Orient))
                .findFirst().orElse(null);
    }

    private static boolean areFacingBreakends(
            final Breakend breakend, final Breakend nextBreakend)
    {
        if(breakend.otherBreakend() == nextBreakend) // ignore DUPs
            return false;

        // order is maintained ie the first must face the second
        if(breakend.Orient == nextBreakend.Orient || !breakend.Chromosome.equals(nextBreakend.Chromosome))
            return false;

        if(breakend.Orient.isReverse())
            return breakend.Position < nextBreakend.Position && nextBreakend.Position - breakend.Position <= PHASED_ASSEMBLY_MAX_TI;
        else
            return breakend.Position > nextBreakend.Position && breakend.Position - nextBreakend.Position <= PHASED_ASSEMBLY_MAX_TI;
    }

    protected static Orientation segmentOrientation(final AlignData alignment, boolean linksEnd)
    {
        return segmentOrientation(alignment.orientation(), linksEnd);
    }

    protected static Orientation segmentOrientation(final Orientation segmentOrientation, boolean linksEnd)
    {
        return (linksEnd == (segmentOrientation.isForward())) ? FORWARD : REVERSE;
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

        relatedAlignments.forEach(x -> altAlignments.addAll(x.allAlignments()));

        return altAlignments;
    }
}
