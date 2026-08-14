package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.tars.common.BwaScoring.maxScoringPrefix;
import static com.hartwig.hmftools.tars.common.BwaScoring.score;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_SCORE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.common.TarsCigarUtils;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

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

    private byte[] refBases(final String chromosome, final int posStart, final int posEnd)
    {
        return BwaScoring.refWindow(mRefGenome, chromosome, posStart, posEnd);
    }

    public boolean enabled()
    {
        return mRefGenome != null;
    }

    public void gateCandidates(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(!enabled() || alignments.isEmpty())
        {
            return;
        }

        ReadBases readBases = ReadBases.of(record);
        if(readBases == null)
        {
            return;
        }

        LiftedAlignment self = alignments.get(0);

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alignment = alignments.get(i);
            if(alignment.LiftedCigar == null)
            {
                continue;
            }

            Outcome outcome = gate(alignment, readBases.forAlignment(alignment));
            if(outcome == null)
            {
                continue;
            }

            if(alignment != self && outcome.dropped())
            {
                alignment.Dropped = true;
                continue;
            }

            if(outcome.alignment() != alignment)
            {
                alignments.set(i, outcome.alignment());
            }
        }
    }

    private record Result(int pos, String cigar, boolean dropped)
    {
    }

    public record Outcome(LiftedAlignment alignment, boolean dropped)
    {
    }

    public Outcome gate(final LiftedAlignment alignment, final byte[] readBases)
    {
        if(alignment == null || alignment.LiftedCigar == null)
        {
            return null;
        }

        Result collapsed = collapseJunctions(
                alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, readBases);
        if(collapsed.pos() == alignment.LiftedPos && collapsed.cigar().equals(alignment.LiftedCigar))
        {
            return new Outcome(alignment, collapsed.dropped());
        }
        return new Outcome(alignment.withLiftedCigar(collapsed.pos(), collapsed.cigar()), collapsed.dropped());
    }

    private Result collapseJunctions(final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefGenome == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        int pos = alignmentStart;
        String working = cigar;

        for(int side = 0; side <= 1; ++side)
        {
            boolean leading = side == 1;

            while(working.indexOf('N') >= 0)
            {
                CollapseResult collapse = tryCollapse(leading, chromosome, pos, working, readBases);
                if(!collapse.collapsed())
                {
                    break;
                }
                pos = collapse.pos();
                working = collapse.cigar();
            }
        }

        boolean dropped = cigar.indexOf('N') >= 0 && working.indexOf('N') < 0;
        return new Result(pos, working, dropped);
    }

    private CollapseResult tryCollapse(
            final boolean leading, final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        return leading
                ? tryCollapseLeading(chromosome, alignmentStart, cigar, readBases)
                : tryCollapseTrailing(chromosome, alignmentStart, cigar, readBases);
    }

    private record Overhang(
            int softclipLength, int anchorLength, int intronIndex, int intronLength, int nearIndex, int window)
    {
    }

    private static Overhang findOverhang(final boolean leading, final String cigar, final List<CigarElement> elements)
    {
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
        {
            return null;
        }

        int outerIndex = leading ? 0 : elements.size() - 1;
        int step = leading ? 1 : -1;

        boolean hasSoftclip = elements.get(outerIndex).getOperator() == CigarOperator.S;
        int softclipLength = hasSoftclip ? elements.get(outerIndex).getLength() : 0;
        int anchorIndex = hasSoftclip ? outerIndex + step : outerIndex;
        int intronIndex = anchorIndex + step;
        int nearIndex = intronIndex + step;
        if(nearIndex < 0 || nearIndex >= elements.size())
        {
            return null;
        }

        CigarElement anchor = elements.get(anchorIndex);
        if(!TarsCigarUtils.isMatchOrEqualOp(anchor.getOperator())
                || !TarsCigarUtils.isMatchOrEqualOp(elements.get(nearIndex).getOperator())
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
        byte[] windowRef = refBases(chromosome, nearEnd + 1, nearEnd + window);
        if(windowRef == null)
        {
            return CollapseResult.unchanged();
        }

        int windowStart = readBases.length - window;
        byte[] anchorRead = Arrays.copyOfRange(readBases, windowStart, windowStart + overhang.anchorLength());
        byte[] anchorRef = refBases(
                chromosome, nearEnd + overhang.intronLength() + 1, nearEnd + overhang.intronLength() + overhang.anchorLength());
        int overhangScore = anchorRef != null ? score(anchorRead, anchorRef) : Integer.MIN_VALUE;
        if(keepJunction(overhangScore, anchorRead, windowRef, overhang.softclipLength()))
        {
            return CollapseResult.unchanged();
        }

        int extended = maxScoringPrefix(Arrays.copyOfRange(readBases, windowStart, readBases.length), windowRef);

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
        byte[] windowRef = refBases(chromosome, nearStart - window, nearStart - 1);
        if(windowRef == null)
        {
            return CollapseResult.unchanged();
        }

        byte[] anchorRead = Arrays.copyOfRange(
                readBases, overhang.softclipLength(), overhang.softclipLength() + overhang.anchorLength());
        byte[] anchorRef = refBases(chromosome, alignmentStart, alignmentStart + overhang.anchorLength() - 1);
        int overhangScore = anchorRef != null ? score(anchorRead, anchorRef) : Integer.MIN_VALUE;
        if(keepJunction(overhangScore, anchorRead, windowRef, overhang.softclipLength()))
        {
            return CollapseResult.unchanged();
        }

        int extended = maxScoringPrefix(reverseArray(Arrays.copyOfRange(readBases, 0, window)), reverseArray(windowRef));

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
        int altRefScore = score(anchorRead, contiguousRef);
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
