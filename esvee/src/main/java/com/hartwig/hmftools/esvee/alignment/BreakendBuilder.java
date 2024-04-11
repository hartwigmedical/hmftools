package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.flipOrientation;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MAX_ZERO_QUALS;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class BreakendBuilder
{
    private final RefGenomeInterface mRefGenome;
    private final AssemblyAlignment mAssemblyAlignment;

    public BreakendBuilder(final RefGenomeInterface refGenome, final AssemblyAlignment assemblyAlignment)
    {
        mRefGenome = refGenome;
        mAssemblyAlignment = assemblyAlignment;
    }

    public void formBreakends(final List<AlignData> alignments, final String fullSequence)
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

        validAlignments.forEach(x -> x.setAdjustedSequenceCoords(fullSequence.length()));

        if(validAlignments.size() == 1)
        {
            if(validAlignments.get(0).maxSoftClipLength() < ALIGNMENT_MIN_SOFT_CLIP)
                return;

            formSingleBreakend(validAlignments.get(0), fullSequence);
        }
        else
        {
            // otherwise handle multiple alignments
            formMultipleBreakends(validAlignments, fullSequence);
        }

        setAlternativeAlignments(zeroQualAlignments);
    }

    private void formSingleBreakend(final AlignData alignment, final String fullSequence)
    {
        int fullSequenceLength = fullSequence.length();

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

        Breakend breakend = new Breakend(alignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, null   );

        breakend.setAnchorLength(alignment.alignedBases());

        mAssemblyAlignment.addBreakend(breakend);
    }

    private void checkOuterSingle(final AlignData alignment, boolean checkStart, final String fullSequence)
    {
        byte orientation = alignment.orientation();

        int softClipLength = checkStart ? alignment.leftSoftClipLength() : alignment.rightSoftClipLength();

        if(softClipLength < ALIGNMENT_MIN_SOFT_CLIP)
            return;

        byte sglOrientation = flipOrientation(orientation);
        int sglPosition = alignment.RefLocation.start();

        int fullSequenceLength = fullSequence.length();
        String insertSequence = checkStart ?
                fullSequence.substring(0, softClipLength) : fullSequence.substring(fullSequenceLength - softClipLength);

        Breakend sglBreakend = new Breakend(
                alignment.RefLocation.Chromosome, sglPosition, sglOrientation, insertSequence, null);

        sglBreakend.setAnchorLength(alignment.alignedBases());

        mAssemblyAlignment.addBreakend(sglBreakend);
    }

    private void formMultipleBreakends(final List<AlignData> alignments, final String fullSequence)
    {
        // look for consecutive non-zero-MQ alignments from which to form breakends
        Collections.sort(alignments, Comparator.comparingInt(x -> x.adjustedSequenceStart()));

        // check for a single at the start or end
        checkOuterSingle(alignments.get(0), true, fullSequence);

        for(int i = 0; i < alignments.size() - 1; ++i)
        {
            AlignData alignment = alignments.get(i);

            byte orientation = alignment.orientation();
            int breakendPosition = alignment.RefLocation.end();

            HomologyData homology = null;
            String insertedBases = "";

            AlignData nextAlignment = alignments.get(i + 1);

            if(alignment.adjustedSequenceEnd() >= nextAlignment.adjustedSequenceStart())
            {
                homology = HomologyData.determineHomology(alignment, nextAlignment, mRefGenome);
            }
            else
            {
                int insertSeqStart = alignment.adjustedSequenceEnd() + 1;
                int insertSeqEnd = nextAlignment.adjustedSequenceStart();
                insertedBases = fullSequence.substring(insertSeqStart, insertSeqEnd);
            }

            Breakend breakend = new Breakend(
                    alignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, homology);

            breakend.setAnchorLength(alignment.alignedBases());

            mAssemblyAlignment.addBreakend(breakend);

            Breakend nextBreakend = new Breakend(
                    nextAlignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, homology);

            nextBreakend.setAnchorLength(nextAlignment.alignedBases());

            mAssemblyAlignment.addBreakend(nextBreakend);
        }

        checkOuterSingle(alignments.get(alignments.size() - 1), false, fullSequence);
    }

    private void setAlternativeAlignments(final List<AlignData> zeroQualAlignments)
    {
        if(zeroQualAlignments.isEmpty() && zeroQualAlignments.size() <= ALIGNMENT_MAX_ZERO_QUALS)
            mAssemblyAlignment.breakends().forEach(x -> x.setAlternativeAlignments(zeroQualAlignments));
    }


    /* OLD LOGIC just for reference sake

                if(alignments.isEmpty())
            {
                assemblyAlignment.setOutcome(NO_RESULT);
                return;
            }

            if(isSequenceLengthMatch(alignments.get(0), fullSequence.length()))
            {
                assemblyAlignment.setOutcome(NON_SV_MATCH);
                return;
            }

            AlignmentOutcome topOutcomeFirst = NO_RESULT;
            AlignmentOutcome topOutcomeSecond = NO_RESULT;

            JunctionAssembly firstAssembly = assemblyAlignment.first();
            JunctionAssembly secondAssembly = assemblyAlignment.second();

            for(AlignData alignment : alignments)
            {
                topOutcomeFirst = assessTopAlignmentResult(firstAssembly, alignment, topOutcomeFirst);

                if(secondAssembly != null)
                    topOutcomeSecond = assessTopAlignmentResult(secondAssembly, alignment, topOutcomeSecond);
            }

            firstAssembly.setAlignmentOutcome(topOutcomeFirst);

            if(secondAssembly != null)
                secondAssembly.setAlignmentOutcome(topOutcomeSecond);


        private static final int ALIGN_MATCH_JUNCTION_POS_BUFFER = 3;

        private static final int ALIGN_MATCH_DIFF_ABS = 2;
        private static final double ALIGN_MATCH_DIFF_PERC = 0.02;
        private static final double ALIGN_MATCH_PERC = 1 - ALIGN_MATCH_DIFF_PERC;

        private boolean isSequenceLengthMatch(final AlignData alignment, final int sequenceLength)
        {
            int fullSequenceLength = sequenceLength;
            int alignmentLength = alignment.RefLocation.baseLength();

            int fullMatchDiff = max(fullSequenceLength - alignmentLength, 0);

            return fullMatchDiff <= ALIGN_MATCH_DIFF_ABS || fullMatchDiff / (double) fullSequenceLength <= ALIGN_MATCH_DIFF_PERC;
        }

        private AlignmentOutcome assessTopAlignmentResult(
                final JunctionAssembly assembly, final AlignData alignment, final AlignmentOutcome existingOutcome)
        {
            AlignmentOutcome outcome = assessAlignmentResult(assembly, alignment);

            if(outcome.exactMatch() && existingOutcome.exactMatch())
                return MULTIPLE;

            return outcome.ordinal() < existingOutcome.ordinal() ? outcome : existingOutcome;
        }

        private AlignmentOutcome assessAlignmentResult(final JunctionAssembly assembly, final AlignData alignment)
        {
            int assemblyRefLength = assembly.refBaseLength();

            boolean sequenceLengthMatch = isSequenceLengthMatch(alignment, assemblyRefLength);

            if(!assembly.junction().chromosome().equals(alignment.RefLocation.Chromosome))
                return sequenceLengthMatch ? ALT_LOC_MATCH : NO_MATCH;

            int assemblyPosStart = assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.junction().Position;
            int assemblyPosEnd = assembly.isForwardJunction() ? assembly.junction().Position : assembly.maxAlignedPosition();

            if(!positionsOverlap(alignment.RefLocation.start(), alignment.RefLocation.end(), assemblyPosStart, assemblyPosEnd))
                return sequenceLengthMatch ? ALT_LOC_MATCH : NO_MATCH;

            boolean matchesJunctionPosition = false;

            if(assembly.isForwardJunction())
            {
                if(abs(alignment.RefLocation.end() - assembly.junction().Position) <= ALIGN_MATCH_JUNCTION_POS_BUFFER)
                    matchesJunctionPosition = true;
            }
            else
            {
                if(abs(alignment.RefLocation.start() - assembly.junction().Position) <= ALIGN_MATCH_JUNCTION_POS_BUFFER)
                    matchesJunctionPosition = true;
            }

            if(!matchesJunctionPosition)
                return sequenceLengthMatch ? ALT_LOC_MATCH : NO_MATCH;

            int refBaseOverlap = min(alignment.RefLocation.end(), assemblyPosEnd) - max(alignment.RefLocation.start(), assemblyPosStart);
            double refBaseOverlapPerc = refBaseOverlap / (double)assemblyRefLength;

            return refBaseOverlapPerc >= ALIGN_MATCH_PERC ? MATCH : PARTIAL;
        }
    }


     */
}
