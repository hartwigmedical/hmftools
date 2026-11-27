package com.hartwig.hmftools.tars.liftback.overhang;

import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_SCORE;
import static com.hartwig.hmftools.tars.liftback.GenomicAlignmentScorer.bestScoringPrefixLength;
import static com.hartwig.hmftools.tars.liftback.GenomicAlignmentScorer.referenceBases;
import static com.hartwig.hmftools.tars.liftback.GenomicAlignmentScorer.ungappedScore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.ReadBases;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class OverhangGate
{
    private final RefGenomeInterface mRefGenome;

    public OverhangGate(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public void gatePlacements(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        ReadBases readBases = new ReadBases(record.getReadBases());
        LiftedAlignment self = alignments.get(0);

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alignment = alignments.get(i);
            byte[] alignmentBases = alignment.forwardStrand() != record.getReadNegativeStrandFlag()
                    ? readBases.getBases() : readBases.getReverseComplement();
            Outcome outcome = gate(alignment, alignmentBases);

            if(alignment != self && outcome.dropped())
            {
                alignments.set(i, alignment.withDropped());
                continue;
            }

            if(outcome.alignment() != alignment)
            {
                alignments.set(i, outcome.alignment());
            }
        }
    }

    public record Outcome(LiftedAlignment alignment, boolean dropped)
    {
    }

    public Outcome gate(final LiftedAlignment alignment, final byte[] readBases)
    {
        int pos = alignment.liftedPosition().Position;
        String working = alignment.liftedCigar();

        for(int side = 0; side <= 1; ++side)
        {
            boolean leading = side == 1;

            while(working.indexOf('N') >= 0)
            {
                CollapseResult collapse = leading
                        ? tryCollapseLeading(alignment.liftedPosition().Chromosome, pos, working, readBases)
                        : tryCollapseTrailing(alignment.liftedPosition().Chromosome, pos, working, readBases);
                if(!collapse.collapsed())
                {
                    break;
                }
                pos = collapse.pos();
                working = collapse.cigar();
            }
        }

        boolean dropped = alignment.cigarHasN() && working.indexOf('N') < 0;
        if(pos == alignment.liftedPosition().Position && working.equals(alignment.liftedCigar()))
        {
            return new Outcome(alignment, dropped);
        }
        return new Outcome(alignment.withLiftedCigar(pos, working), dropped);
    }

    private record Overhang(
            int softclipLength, int anchorLength, int intronIndex, int intronLength, int nearIndex, int window)
    {
    }

    private static Overhang findOverhang(final boolean leading, final String cigar, final List<CigarElement> elements)
    {
        int outerIndex = leading ? 0 : elements.size() - 1;
        int step = leading ? 1 : -1;

        boolean hasSoftclip = elements.get(outerIndex).getOperator() == CigarOperator.S;
        int softclipLength = hasSoftclip ? elements.get(outerIndex).getLength() : 0;
        int anchorIndex = hasSoftclip ? outerIndex + step : outerIndex;
        int intronIndex = anchorIndex + step;
        int nearIndex = intronIndex + step;
        CigarElement anchor = elements.get(anchorIndex);
        if(!CigarUtils.isMatchOrEqualOp(anchor.getOperator())
                || !CigarUtils.isMatchOrEqualOp(elements.get(nearIndex).getOperator())
                || elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return null;
        }

        if(softclipLength == 0 && cigar.indexOf('N') == cigar.lastIndexOf('N'))
        {
            return null;
        }

        if(anchor.getLength() > MIN_OVERHANG_LENGTH)
        {
            return null;
        }

        return new Overhang(
                softclipLength, anchor.getLength(), intronIndex, elements.get(intronIndex).getLength(), nearIndex,
                anchor.getLength() + softclipLength);
    }

    private CollapseResult tryCollapseTrailing(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        Overhang overhang = findOverhang(false, cigar, elements);
        if(overhang == null)
        {
            return CollapseResult.unchanged();
        }

        int window = overhang.window();
        int nearEnd = alignmentStart + CigarUtils.cigarAlignedLength(elements.subList(0, overhang.intronIndex())) - 1;
        byte[] windowRef = referenceBases(mRefGenome, chromosome, nearEnd + 1, nearEnd + window);

        int windowStart = readBases.length - window;
        byte[] anchorRead = Arrays.copyOfRange(readBases, windowStart, windowStart + overhang.anchorLength());
        byte[] anchorRef = referenceBases(
                mRefGenome, chromosome, nearEnd + overhang.intronLength() + 1,
                nearEnd + overhang.intronLength() + overhang.anchorLength());
        int overhangScore = ungappedScore(anchorRead, anchorRef);
        if(keepJunction(overhangScore, anchorRead, windowRef, overhang.softclipLength()))
        {
            return CollapseResult.unchanged();
        }

        int extended = bestScoringPrefixLength(Arrays.copyOfRange(readBases, windowStart, readBases.length), windowRef);

        CigarElement nearExon = elements.get(overhang.nearIndex());
        List<CigarElement> merged = new ArrayList<>(elements.subList(0, overhang.nearIndex()));
        merged.add(new CigarElement(nearExon.getLength() + extended, CigarOperator.M));
        if(window > extended)
        {
            merged.add(new CigarElement(window - extended, CigarOperator.S));
        }
        return new CollapseResult(true, alignmentStart, CigarUtils.cigarElementsToStr(merged));
    }

    private CollapseResult tryCollapseLeading(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        Overhang overhang = findOverhang(true, cigar, elements);
        if(overhang == null)
        {
            return CollapseResult.unchanged();
        }

        int window = overhang.window();
        int nearStart = alignmentStart + overhang.anchorLength() + overhang.intronLength();
        byte[] windowRef = referenceBases(mRefGenome, chromosome, nearStart - window, nearStart - 1);

        byte[] anchorRead = Arrays.copyOfRange(
                readBases, overhang.softclipLength(), overhang.softclipLength() + overhang.anchorLength());
        byte[] anchorRef = referenceBases(mRefGenome, chromosome, alignmentStart, alignmentStart + overhang.anchorLength() - 1);
        int overhangScore = ungappedScore(anchorRead, anchorRef);
        if(keepJunction(overhangScore, anchorRead, windowRef, overhang.softclipLength()))
        {
            return CollapseResult.unchanged();
        }

        int extended = bestScoringPrefixLength(reverseArray(Arrays.copyOfRange(readBases, 0, window)), reverseArray(windowRef));

        CigarElement nearExon = elements.get(overhang.nearIndex());
        List<CigarElement> merged = new ArrayList<>();
        if(window > extended)
        {
            merged.add(new CigarElement(window - extended, CigarOperator.S));
        }
        merged.add(new CigarElement(extended + nearExon.getLength(), CigarOperator.M));
        merged.addAll(elements.subList(overhang.nearIndex() + 1, elements.size()));
        return new CollapseResult(true, nearStart - extended, CigarUtils.cigarElementsToStr(merged));
    }

    private static boolean keepJunction(final int overhangScore, final byte[] anchorRead, final byte[] windowRef, final int softclipLength)
    {
        if(softclipLength > 0)
        {
            return overhangScore > MIN_OVERHANG_SCORE;
        }
        if(overhangScore > 0)
        {
            return true;
        }
        byte[] contiguousRef = Arrays.copyOfRange(windowRef, 0, anchorRead.length);
        int altRefScore = ungappedScore(anchorRead, contiguousRef);
        return altRefScore <= overhangScore;
    }

    private record CollapseResult(boolean collapsed, int pos, String cigar)
    {
        static CollapseResult unchanged()
        {
            return new CollapseResult(false, 0, null);
        }
    }

}
