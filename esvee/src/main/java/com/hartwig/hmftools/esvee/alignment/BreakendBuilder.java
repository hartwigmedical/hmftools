package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
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
            if(validAlignments.get(0).maxSoftClipLength() < ALIGNMENT_MIN_SOFT_CLIP)
                return;

            formSingleBreakend(validAlignments.get(0));
        }
        else
        {
            // otherwise handle multiple alignments
            formMultipleBreakends(validAlignments);
        }

        setAlternativeAlignments(zeroQualAlignments);
    }

    private void formSingleBreakend(final AlignData alignment)
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

        Breakend breakend = new Breakend(alignment.RefLocation.Chromosome, breakendPosition, orientation, insertedBases, null);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), fullSequenceLength, alignment.sequenceStart(), orientation, 0, alignment);

        breakend.addSegment(segment);

        mAssemblyAlignment.addBreakend(breakend);
    }

    private boolean checkOuterSingle(final AlignData alignment, boolean checkStart, int nextSegmentIndex)
    {
        byte sglOrientation = segmentOrientation(alignment, checkStart);

        int softClipLength = checkStart ? alignment.leftSoftClipLength() : alignment.rightSoftClipLength();

        if(softClipLength < ALIGNMENT_MIN_SOFT_CLIP)
            return false;

        String fullSequence = mAssemblyAlignment.fullSequence();
        int fullSequenceLength = mAssemblyAlignment.fullSequenceLength();

        int sglPosition = checkStart ? alignment.RefLocation.start() : alignment.RefLocation.end();

        String insertSequence = checkStart ?
                fullSequence.substring(0, softClipLength) : fullSequence.substring(fullSequenceLength - softClipLength);

        Breakend breakend = new Breakend(
                alignment.RefLocation.Chromosome, sglPosition, sglOrientation, insertSequence, null);

        BreakendSegment segment = new BreakendSegment(
                mAssemblyAlignment.id(), fullSequenceLength, alignment.sequenceStart(), sglOrientation, nextSegmentIndex, alignment);

        breakend.addSegment(segment);

        mAssemblyAlignment.addBreakend(breakend);

        return true;
    }

    private void formMultipleBreakends(final List<AlignData> alignments)
    {
        // look for consecutive non-zero-MQ alignments from which to form breakends
        Collections.sort(alignments, Comparator.comparingInt(x -> x.sequenceStart()));

        String fullSequence = mAssemblyAlignment.fullSequence();

        int nextSegmentIndex = 0;

        // check for a single at the start or end
        if(checkOuterSingle(alignments.get(0), true, nextSegmentIndex))
            nextSegmentIndex++;

        for(int i = 0; i < alignments.size() - 1; ++i)
        {
            AlignData alignment = alignments.get(i);

            byte breakendOrientation = segmentOrientation(alignment, true);
            int breakendPosition = alignment.RefLocation.end();

            HomologyData homology = null;
            String insertedBases = "";

            AlignData nextAlignment = alignments.get(i + 1);

            if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
            {
                homology = HomologyData.determineHomology(alignment, nextAlignment, mRefGenome);
            }
            else
            {
                int insertSeqStart = alignment.sequenceEnd() + 1;
                int insertSeqEnd = nextAlignment.sequenceStart();
                insertedBases = fullSequence.substring(insertSeqStart, insertSeqEnd);
            }

            Breakend breakend = new Breakend(
                    alignment.RefLocation.Chromosome, breakendPosition, breakendOrientation, insertedBases, homology);

            mAssemblyAlignment.addBreakend(breakend);

            BreakendSegment segment = new BreakendSegment(
                    mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), alignment.sequenceEnd(),
                    breakendOrientation, nextSegmentIndex++, alignment);

            breakend.addSegment(segment);

            byte nextOrientation = segmentOrientation(nextAlignment, false);
            int nextPosition = nextAlignment.RefLocation.start();

            Breakend nextBreakend = new Breakend(
                    nextAlignment.RefLocation.Chromosome, nextPosition, nextOrientation, insertedBases, homology);

            mAssemblyAlignment.addBreakend(nextBreakend);

            BreakendSegment nextSegment = new BreakendSegment(
                    mAssemblyAlignment.id(), mAssemblyAlignment.fullSequenceLength(), nextAlignment.sequenceStart(),
                    nextOrientation, nextSegmentIndex++, nextAlignment);

            nextBreakend.addSegment(nextSegment);

            breakend.setOtherBreakend(nextBreakend);
            nextBreakend.setOtherBreakend(breakend);
        }

        checkOuterSingle(alignments.get(alignments.size() - 1), false, nextSegmentIndex);
    }

    protected static byte segmentOrientation(final AlignData alignment, boolean linksEnd)
    {
        return (linksEnd == (alignment.orientation() == POS_ORIENT)) ? POS_ORIENT : NEG_ORIENT;
    }

    private void setAlternativeAlignments(final List<AlignData> zeroQualAlignments)
    {
        if(zeroQualAlignments.isEmpty() && zeroQualAlignments.size() <= ALIGNMENT_MAX_ZERO_QUALS)
            mAssemblyAlignment.breakends().forEach(x -> x.setAlternativeAlignments(zeroQualAlignments));
    }
}
