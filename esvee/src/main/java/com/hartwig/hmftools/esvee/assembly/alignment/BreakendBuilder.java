package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineSequenceCount;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.common.saga.SagaAlignment;
import com.hartwig.hmftools.esvee.common.saga.SagaAssembly;
import com.hartwig.hmftools.esvee.common.saga.SagaBreakend;
import com.hartwig.hmftools.esvee.common.saga.SagaVariant;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

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

    public boolean formBreakendsFromSaga(final SagaAlignment sagaAlignment, final AlignData alignData)
    {
        // For a SAGA-matched variant, use the breakends stated in the SAGA resource. Why?
        //   - Some cases where homology causes alignment to be misleading, resulting in wrong breakends.
        //   - Recover SGLs to full variants given that one half matches a SAGA variant.

        // Only use SAGA breakends for the simple case of 0 or 1 links. Chained assemblies require more care, not implemented for now.
        boolean isChained = mAssemblyAlignment.phaseSet() != null && mAssemblyAlignment.phaseSet().assemblyLinks().size() > 1;
        if(isChained)
        {
            return false;
        }

        // Reverse strand match shouldn't occur and isn't handled here.
        if(!sagaAlignment.isForward())
        {
            return false;
        }

        SagaAssembly sagaAssembly = sagaAlignment.sagaAssembly();
        SagaVariant sagaVariant = sagaAssembly.variant();
        SagaBreakend lowerSagaBreakend = sagaVariant.breakend1();
        SagaBreakend upperSagaBreakend = sagaVariant.breakend2();

        // SAGA variants should only be indels, but check because other types are not handled here.
        if(!lowerSagaBreakend.chromosome().equals(upperSagaBreakend.chromosome()))
        {
            return false;
        }

        String rawInsertSequence;
        int indelLength;
        boolean lowerBreakendInferred;
        boolean upperBreakendInferred;
        if(sagaVariant.isInsert())
        {
            // Insertion - the assembly overlaps the SAGA insert sequence, which is bounded by the junction offsets.

            // Replace the ESVEE-assembled insert sequence with the SAGA insert sequence. This is because:
            //   - In the case of 2 unlinked SGLs recovered to paired breakends, deduping requires the insert sequence to be the same.
            //   - The differences between the phased assembly and SAGA sequence are known to be small, otherwise no SAGA match.
            //   - Allows more advanced resolution of the SAGA breakends (e.g. template insertions) to be precomputed, in the future.
            rawInsertSequence = sagaVariant.insertSequence();

            indelLength = rawInsertSequence.length();

            // If the phased assembly didn't span the whole SAGA variant, take note that the other breakend is inferred from SAGA.
            List<Integer> sagaJunctionOffsets = sagaAssembly.junctionOffsets();
            int startInferredLength = sagaAlignment.sagaStart() - sagaJunctionOffsets.get(0);
            int endInferredLength = sagaJunctionOffsets.get(1) - sagaAlignment.sagaEnd();
            lowerBreakendInferred = startInferredLength > 0;
            upperBreakendInferred = endInferredLength > 0;
        }
        else
        {
            // Deletion - no insert sequence and the assembly spans the junction (otherwise it couldn't have matched the SAGA variant).

            rawInsertSequence = "";
            // Using the SAGA breakends can differ in a few bases of deletion, which is ok. Similar explanation to the insert case above.
            indelLength = upperSagaBreakend.position() - lowerSagaBreakend.position() - 1;
            lowerBreakendInferred = false;
            upperBreakendInferred = false;
        }

        // Calculate the indices of the (adjusted) breakends within the phased assembly sequence.
        // Note these can be out of bounds if the assembly doesn't span the whole SAGA variant.
        // These are only used for fragment support counting.
        List<Integer> assemblyJunctionOffsets = sagaAlignment.queryJunctionOffsets();
        int assemblyIndelStart = assemblyJunctionOffsets.get(0);
        int assemblyIndelEnd = assemblyJunctionOffsets.size() > 1 ? assemblyJunctionOffsets.get(1) : assemblyIndelStart;

        IndelCoords indelCoords = new IndelCoords(lowerSagaBreakend.position(), upperSagaBreakend.position(), indelLength);
        indelCoords.setInsertedBases(rawInsertSequence);

        IndelBreakendData breakendData =
                calcIndelBreakendsWithHomology(lowerSagaBreakend.chromosome(), indelCoords, assemblyIndelStart, assemblyIndelEnd);

        Pair<Breakend, Breakend> breakends = createIndelBreakends(alignData, breakendData);
        Breakend lowerBreakend = breakends.getLeft();
        Breakend upperBreakend = breakends.getRight();

        if(lowerBreakendInferred)
        {
            lowerBreakend.setSagaInferred();
        }
        if(upperBreakendInferred)
        {
            upperBreakend.setSagaInferred();
        }

        return true;
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

        setMaxLocalRepeat();
    }

    private boolean formIndelBreakends(final AlignData alignment)
    {
        IndelCoords indelCoords = getIndelAlignmentIndelCoords(alignment);

        if(indelCoords == null || indelCoords.Length < MIN_INDEL_LENGTH)
            return false;

        int assemblyIndelStart = alignment.sequenceStart()
                + getReadIndexFromPosition(alignment.positionStart(), alignment.cigarElements(), indelCoords.PosStart, true, false)
                // Advance from the last ref base to the first indel base.
                + 1
                // back out soft-clip since the method above includes that
                - leftSoftClipLength(alignment.cigarElements());
        int assemblyIndelEnd = assemblyIndelStart + (indelCoords.isInsert() ? indelCoords.Length : 0);

        boolean isChained = mAssemblyAlignment.phaseSet() != null && mAssemblyAlignment.phaseSet().assemblyLinks().size() > 1;
        if(!isChained)
        {
            // the indel must have sufficient bases either side of it to be called
            int leftAnchorLength = assemblyIndelStart;
            int rightAnchorLength = alignment.segmentLength() - assemblyIndelEnd;

            if(leftAnchorLength < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH || rightAnchorLength < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH)
                return false;
        }

        String assemblySeq = mAssemblyAlignment.fullSequence();
        if(indelCoords.isInsert() && assemblyIndelStart >= 0 && assemblyIndelEnd <= assemblySeq.length())
        {
            indelCoords.setInsertedBases(assemblySeq.substring(assemblyIndelStart, assemblyIndelEnd));
        }

        IndelBreakendData breakendData =
                calcIndelBreakendsWithHomology(alignment.chromosome(), indelCoords, assemblyIndelStart, assemblyIndelEnd);

        createIndelBreakends(alignment, breakendData);

        return true;
    }

    @Nullable
    private IndelCoords getIndelAlignmentIndelCoords(final AlignData alignment)
    {
        // parse the CIGAR to get the indel coords, first looking for a close match to what was expected
        IndelCoords indelCoords = null;
        if(mAssemblyAlignment.phaseSet() != null)
        {
            if(mAssemblyAlignment.phaseSet().assemblyLinks().size() == 1)
            {
                AssemblyLink assemblyLink = mAssemblyAlignment.phaseSet().assemblyLinks().get(0);
                StructuralVariantType svType = assemblyLink.svType();

                if(svType == DEL || svType == DUP || svType == INS)
                {
                    int svLength = assemblyLink.length();

                    if(svLength == 0 && !assemblyLink.insertedBases().isEmpty()) // same-base DUP/INS, change implied Cigar insert length
                    {
                        svLength = assemblyLink.insertedBases().length();
                        svType = INS;
                    }

                    CigarElement specificIndel = new CigarElement(svLength, svType == DEL ? D : I);
                    indelCoords = IndelCoords.findMatchingIndelCoords(alignment.positionStart(), alignment.cigarElements(), specificIndel);
                }
            }
        }

        // if no match was found, just take the longest
        if(indelCoords == null)
        {
            indelCoords = findIndelCoords(alignment.positionStart(), alignment.cigarElements(), MIN_INDEL_LENGTH);
        }

        return indelCoords;
    }

    private record IndelBreakendData(
            String chromosome,
            int lowerPosition,
            Orientation lowerOrient,
            int upperPosition,
            Orientation upperOrient,
            String insertedBases,
            @Nullable
            HomologyData homology,
            int assemblyIndelStart,
            int assemblyIndelEnd
    )
    {
    }

    private IndelBreakendData calcIndelBreakendsWithHomology(final String chromosome, final IndelCoords indelCoords, int assemblyIndelStart,
            int assemblyIndelEnd)
    {
        int lowerPosition = indelCoords.PosStart;
        int upperPosition = indelCoords.PosEnd;

        Orientation lowerOrient = FORWARD;
        Orientation upperOrient = REVERSE;

        String insertedBases = indelCoords.insertedBases();

        HomologyData homology = null;

        if(indelCoords.isInsert())
        {
            if(!insertedBases.isEmpty())
            {
                String basesStart = insertedBases;

                int homPosEnd = upperPosition;
                String basesEnd = mRefGenome.getBaseString(chromosome, homPosEnd, homPosEnd + insertedBases.length() - 1);

                homology = HomologyData.determineIndelHomology(basesStart, basesEnd, insertedBases.length());
                assert homology != null;

                if(!homology.Homology.isEmpty())
                {
                    // convert an INS to a DUP and reassess homology
                    int totalInexactHomology = homology.InexactEnd - homology.InexactStart;

                    assemblyIndelStart += totalInexactHomology;
                    lowerPosition += 1;
                    upperPosition += totalInexactHomology - 1;
                    lowerOrient = REVERSE;
                    upperOrient = FORWARD;
                    int newHomologyLength = 0;
                    if(totalInexactHomology == insertedBases.length())
                    {
                        String postIndelBases =
                                mRefGenome.getBaseString(chromosome, upperPosition + 1, upperPosition + insertedBases.length());

                        while(newHomologyLength < insertedBases.length() && newHomologyLength < postIndelBases.length()
                                && insertedBases.charAt(newHomologyLength) == postIndelBases.charAt(newHomologyLength))
                        {
                            ++newHomologyLength;
                        }
                    }

                    if(newHomologyLength > 0)
                    {
                        String newHomologyBases = insertedBases.substring(0, newHomologyLength);

                        homology = new HomologyData(newHomologyBases, 0, newHomologyLength, 0, newHomologyLength);
                    }
                    else
                    {
                        homology = null;
                    }

                    insertedBases = insertedBases.substring(totalInexactHomology);
                }
            }
        }
        else
        {
            // check for exact homology at the bases either side of the delete
            int maxLength = indelCoords.Length - 1;
            int homPosStart = lowerPosition + 1;
            String basesStart = mRefGenome.getBaseString(chromosome, homPosStart, homPosStart + maxLength);

            int homPosEnd = upperPosition;
            String basesEnd = mRefGenome.getBaseString(chromosome, homPosEnd, homPosEnd + maxLength);

            homology = HomologyData.determineIndelHomology(basesStart, basesEnd, maxLength);
            assert homology != null;

            if(homology.Homology.isEmpty())
            {
                homology = null;
            }
            else
            {
                // Adjust the breakends to be middle-aligned. Doesn't change the nature of the variant, however.
                lowerPosition += abs(homology.ExactStart);
                upperPosition += abs(homology.ExactStart);
            }
        }

        return new IndelBreakendData(
                chromosome, lowerPosition, lowerOrient, upperPosition, upperOrient, insertedBases, homology, assemblyIndelStart, assemblyIndelEnd);
    }

    private Pair<Breakend, Breakend> createIndelBreakends(final AlignData alignData, final IndelBreakendData breakendData)
    {
        Breakend lowerBreakend = new Breakend(
                mAssemblyAlignment, breakendData.chromosome(), breakendData.lowerPosition(), breakendData.lowerOrient(),
                breakendData.insertedBases(), breakendData.homology());

        Breakend upperBreakend = new Breakend(
                mAssemblyAlignment, breakendData.chromosome(), breakendData.upperPosition(), breakendData.upperOrient(),
                breakendData.insertedBases(), breakendData.homology());

        int[] indelSequenceIndices = new int[] { breakendData.assemblyIndelStart(), breakendData.assemblyIndelEnd() };
        BreakendSegment segment =
                new BreakendSegment(mAssemblyAlignment.id(), alignData.sequenceStart(), 0, alignData, indelSequenceIndices);
        lowerBreakend.addSegment(segment);
        upperBreakend.addSegment(segment);

        lowerBreakend.setOtherBreakend(upperBreakend);
        upperBreakend.setOtherBreakend(lowerBreakend);

        mAssemblyAlignment.addBreakend(lowerBreakend);
        mAssemblyAlignment.addBreakend(upperBreakend);

        return Pair.of(lowerBreakend, upperBreakend);
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
            breakendPosition = alignment.positionStart();
            orientation = REVERSE;
            softClipLength = alignment.leftSoftClipLength();
        }
        else
        {
            breakendPosition = alignment.positionEnd();
            orientation = FORWARD;
            softClipLength = alignment.rightSoftClipLength();
        }

        if(alignment.orientation().isForward())
        {
            if(orientation.isReverse())
                insertedBases = fullSequence.substring(0, softClipLength);
            else
                insertedBases = fullSequence.substring(fullSequenceLength - softClipLength);
        }
        else
        {
            if(orientation.isReverse())
                insertedBases = fullSequence.substring(fullSequenceLength - softClipLength);
            else
                insertedBases = fullSequence.substring(0, softClipLength);

            insertedBases = reverseComplementBases(insertedBases);
        }

        int requiredSoftClipLength = calcRequiredSoftClipLength(insertedBases, orientation);

        if(softClipLength < requiredSoftClipLength)
            return;

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.chromosome(), breakendPosition, orientation, insertedBases, null);

        BreakendSegment segment = new BreakendSegment(mAssemblyAlignment.id(), alignment.sequenceStart(), 0, alignment);

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

        int requiredSoftClipLength = calcRequiredSoftClipLength(insertSequence, sglOrientation);

        if(softClipLength < requiredSoftClipLength)
            return false;

        Breakend breakend = new Breakend(
                mAssemblyAlignment, alignment.chromosome(), sglPosition, sglOrientation, insertSequence, null);

        BreakendSegment segment = new BreakendSegment(mAssemblyAlignment.id(), alignment.sequenceStart(), nextSegmentIndex, alignment);

        breakend.addSegment(segment);

        mAssemblyAlignment.addBreakend(breakend);

        List<AlternativeAlignment> altAlignments = checkStart ?
                buildAltAlignments(alignment, null, zeroQualAlignments) : buildAltAlignments(null, alignment, zeroQualAlignments);

        if(altAlignments != null && !altAlignments.isEmpty())
            breakend.setAlternativeAlignments(altAlignments);

        return true;
    }

    private int calcRequiredSoftClipLength(final String insertSequence, final Orientation orientation)
    {
        if(isLineInsertion(insertSequence, orientation))
        {
            return LINE_MIN_EXTENSION_LENGTH;
        }
        else if(mAssemblyAlignment.isSagaMatched())
        {
            return ALIGNMENT_MIN_SOFT_CLIP_LOWER;
        }
        else
        {
            return ALIGNMENT_MIN_SOFT_CLIP;
        }
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
                homology = HomologyData.determineHomology(fullSequence, alignment, nextAlignment);

                if(homology == null)
                {
                    SV_LOGGER.debug("assembly({}) failed to determine homology", mAssemblyAlignment);
                }

                if(alignment.isReverse() && nextAlignment.isReverse())
                {
                    homology = homology.invert(true, false);
                }
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

            String assemblyInsertedBases = breakendOrientation.isForward() ? insertedBases : reverseComplementBases(insertedBases);

            Breakend breakend = new Breakend(
                    mAssemblyAlignment, breakendChr, breakendPosition, breakendOrientation, assemblyInsertedBases, firstHomology);

            mAssemblyAlignment.addBreakend(breakend);

            BreakendSegment segment =
                    new BreakendSegment(mAssemblyAlignment.id(), alignment.sequenceStart(), nextSegmentIndex++, alignment);

            breakend.addSegment(segment);

            String nextInsertedBases = nextOrientation.isReverse() ? insertedBases : reverseComplementBases(insertedBases);

            Breakend nextBreakend = new Breakend(
                    mAssemblyAlignment, nextChr, nextPosition, nextOrientation, nextInsertedBases, nextHomology);

            mAssemblyAlignment.addBreakend(nextBreakend);

            BreakendSegment nextSegment = new BreakendSegment(
                    mAssemblyAlignment.id(), nextAlignment.sequenceStart(), nextSegmentIndex++, nextAlignment);

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
            else if(nextBreakend.svType() == DUP && nextBreakend.Orient.isReverse() && i < mAssemblyAlignment.breakends().size() - 2)
            {
                // handle links to the upper breakend of a short DUP
                Breakend dupOtherBreakend = mAssemblyAlignment.breakends().get(i + 2);

                if(dupOtherBreakend == nextBreakend.otherBreakend() && areFacingBreakends(breakend, dupOtherBreakend))
                {
                    breakend.addFacingBreakend(dupOtherBreakend);
                    dupOtherBreakend.addFacingBreakend(breakend);
                }
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
                    .filter(x -> positionsOverlap(x.sequenceStart(), x.sequenceEnd(), alignmentStart.sequenceStart(), alignmentEnd.sequenceStart()))
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

    private void setMaxLocalRepeat()
    {
        // maximum repeat at the junction is used by the caller for filtering
        if(mAssemblyAlignment.breakends().size() != 2)
            return;

        if(!mAssemblyAlignment.breakends().stream().allMatch(x -> x.isShortLocalDelDupIns()))
            return;

        int maxRepeatLength = 0;

        for(JunctionAssembly assembly : mAssemblyAlignment.assemblies())
        {
            if(assembly.repeatInfo() != null)
            {
                int maxRepeat = assembly.repeatInfo().stream().mapToInt(x -> x.Count).max().orElse(0);
                maxRepeatLength = max(maxRepeatLength, maxRepeat);
            }
        }

        int maxLocalRepeat = maxRepeatLength;
        mAssemblyAlignment.breakends().forEach(x -> x.setMaxLocalRepeat(maxLocalRepeat));
    }
}
