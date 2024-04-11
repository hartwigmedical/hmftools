package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_SOFT_CLIP;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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

    public void formBreakends(final List<AlignData> alignments)
    {
        if(alignments.isEmpty())
            return;

        if(alignments.size() == 1)
        {
            if(alignments.get(0).maxSoftClipLength() < ALIGNMENT_MIN_SOFT_CLIP)
                return;

            formSingleBreakend(alignments);
            return;
        }

        int nonZeroMapQualCount = (int)alignments.stream().filter(x -> x.MapQual == 0).count();

        if(alignments.stream().allMatch(x -> x.MapQual == 0))
            return;

        // otherwise handle multiple alignments

    }

    private void formSingleBreakend(final List<AlignData> alignments)
    {
        AlignData alignment = alignments.get(0);

        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = fullSequence.length();

        int breakendPosition;
        byte orientation;
        int softClipLength;
        String insertedBases;

        if(alignment.SoftClipLeft >= alignment.SoftClipRight)
        {
            breakendPosition = alignment.RefLocation.start();
            orientation = NEG_ORIENT;
            softClipLength = alignment.SoftClipLeft;
            insertedBases = fullSequence.substring(0, softClipLength);
        }
        else
        {
            breakendPosition = alignment.RefLocation.end();
            orientation = POS_ORIENT;
            softClipLength = alignment.SoftClipRight;
            insertedBases = fullSequence.substring(fullSequenceLength - softClipLength);
        }

        Breakend breakend = new Breakend(alignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, SGL);

        mAssemblyAlignment.addBreakend(breakend);
    }

    private void formMultipleBreakends(final List<AlignData> alignments)
    {
        // look for consecutive non-zero-MQ alignments from which to form breakends
        Collections.sort(alignments, Comparator.comparingInt(x -> x.SequenceStart));



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
