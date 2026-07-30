package com.hartwig.hmftools.tars.liftback.overhang;

import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_EXTEND;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_OPEN;
import static com.hartwig.hmftools.tars.common.TarsConstants.MATCH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MISMATCH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_SCORE;

import static htsjdk.samtools.util.SequenceUtil.basesEqual;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Re-scores short terminal splice-junction overhangs against the reference and collapses unsupported ones; README Step 1.
// TODO revisit this later to simplify
public class OverhangGate
{
    private final RefGenomeInterface mRefGenome;

    public OverhangGate(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    // Null unless the whole requested window is present. Every pass here treats a partial window as no reference, and
    // a query overshooting a contig end must not be scored on truncated bases - which is also what the underlying
    // implementations disagree on (RefGenomeSource lets htsjdk throw past the end, MockRefGenome returns empty).
    private byte[] refBases(final String chromosome, final int posStart, final int posEnd)
    {
        if(mRefGenome == null || posStart < 1 || posEnd < posStart || posEnd > mRefGenome.getChromosomeLength(chromosome))
        {
            return null;
        }

        byte[] bases = mRefGenome.getBases(chromosome, posStart, posEnd);
        return bases != null && bases.length == posEnd - posStart + 1 ? bases : null;
    }

    public boolean enabled()
    {
        return mRefGenome != null;
    }

    // Final (start, cigar) of a gate pass, plus whether an XA-alt placement should be dropped from the tag.
    public record Result(int pos, String cigar, boolean dropped)
    {
    }

    public Result collapseJunctions(final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefGenome == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        int pos = alignmentStart;
        String working = cigar;

        // Trailing side first, then leading. Each side loops because collapsing one junction can expose another behind
        // it; a side stops as soon as its outermost junction is kept.
        for(int side = 0; side <= 1; ++side)
        {
            boolean leading = side == 1;

            while(working.indexOf('N') >= 0)
            {
                CollapseResult collapse = leading
                        ? tryCollapseLeading(chromosome, pos, working, readBases)
                        : tryCollapseTrailing(chromosome, pos, working, readBases);
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

    // The gate's shape: an optional outermost softclip, the anchor exon, its intron, and the exon the anchor would fold
    // into. Null when the cigar is not that shape, or when this junction is not a candidate for collapsing at all.
    private record Overhang(
            int softclipLength, int anchorLength, int intronIndex, int intronLength, int nearIndex, int window)
    {
    }

    private static Overhang findOverhang(final boolean leading, final List<CigarElement> elements)
    {
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
        {
            return null;
        }

        // shapes: "<...>M nN yM [eS]" trailing, "[eS] yM nN <...>M" leading
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
        if(!isMatch(anchor.getOperator()) || !isMatch(elements.get(nearIndex).getOperator())
                || elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return null;
        }

        // a lone junction with no adjacent clip is bwa's own call, not a lift artifact, so it is never gated
        if(softclipLength == 0 && countIntrons(elements) <= 1)
        {
            return null;
        }

        // a long overhang is trusted outright (a real junction even with a couple of mismatches); only short ones are scored
        if(anchor.getLength() > MIN_OVERHANG_LENGTH)
        {
            return null;
        }

        return new Overhang(
                softclipLength, anchor.getLength(), intronIndex, elements.get(intronIndex).getLength(), nearIndex,
                anchor.getLength() + softclipLength);
    }

    // Folds a trailing anchor back into the exon before its intron, keeping whatever of the window matches the genome.
    private CollapseResult tryCollapseTrailing(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        Overhang overhang = findOverhang(false, elements);
        if(overhang == null)
        {
            return CollapseResult.unchanged();
        }

        int window = overhang.window();
        // genomic end of the near exon, just before the intron
        int nearEnd = alignmentStart + CigarUtils.cigarAlignedLength(elements.subList(0, overhang.intronIndex())) - 1;
        byte[] windowRef = refBases(chromosome, nearEnd + 1, nearEnd + window);
        if(windowRef == null)
        {
            return CollapseResult.unchanged();
        }

        // score the anchor where the junction would place it, on its own exon
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

    // The mirror of tryCollapseTrailing. The shape check is shared; only the genomic coordinates and the rebuilt cigar
    // differ, and both are clearer written out than hidden behind a direction abstraction.
    private CollapseResult tryCollapseLeading(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        Overhang overhang = findOverhang(true, elements);
        if(overhang == null)
        {
            return CollapseResult.unchanged();
        }

        int window = overhang.window();
        // a leading softclip consumes no reference, so the near exon starts just past anchor + intron
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

        // the near-exon boundary sits at the END of a leading window, so reverse both to walk the score outward
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

    // Whether a short overhang's junction survives the gate.
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

    // Result of one collapse step.
    private record CollapseResult(boolean collapsed, int pos, String cigar)
    {
        static CollapseResult unchanged()
        {
            return new CollapseResult(false, 0, null);
        }
    }

    public Result extendSoftClips(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefGenome == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.isEmpty() || CigarUtils.hasHardClip(elements))
        {
            return new Result(alignmentStart, cigar, false);
        }

        List<CigarElement> working = new ArrayList<>(elements);
        int alignmentEnd = alignmentStart + CigarUtils.cigarAlignedLength(working) - 1;

        int trailExtended = extendSoftClip(false, chromosome, alignmentStart, alignmentEnd, readBases, working);
        int leadExtended = extendSoftClip(true, chromosome, alignmentStart, alignmentEnd, readBases, working);

        if(trailExtended == 0 && leadExtended == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        return new Result(alignmentStart - leadExtended, CigarUtils.cigarElementsToStr(working), false);
    }

    // Walks a terminal softclip outward from its M boundary, converting the bases that match the genome back to M.
    private int extendSoftClip(
            final boolean leading, final String chromosome, final int alignmentStart, final int alignmentEnd,
            final byte[] readBases, final List<CigarElement> working)
    {
        if(working.size() < 2)
        {
            return 0;
        }

        int softclipIndex = leading ? 0 : working.size() - 1;
        int matchedIndex = leading ? 1 : working.size() - 2;

        CigarElement softclip = working.get(softclipIndex);
        if(softclip.getOperator() != CigarOperator.S || !isMatchedOp(working.get(matchedIndex).getOperator()))
        {
            return 0;
        }

        // a leading clip walks back from the alignment start, so it cannot run past the start of the chromosome
        int budget = leading ? Math.min(softclip.getLength(), alignmentStart - 1) : softclip.getLength();
        if(budget <= 0)
        {
            return 0;
        }

        byte[] refBases = leading
                ? refBases(chromosome, alignmentStart - budget, alignmentStart - 1)
                : refBases(chromosome, alignmentEnd + 1, alignmentEnd + budget);
        if(refBases == null || refBases.length == 0)
        {
            return 0;
        }

        int walkLength = Math.min(refBases.length, budget);
        int softclipLength = softclip.getLength();
        int extended;
        if(leading)
        {
            // a leading softclip's M boundary is at its inner end, so reverse both windows to walk boundary-outward
            extended = maxScoringPrefix(
                    reverseArray(Arrays.copyOfRange(readBases, softclipLength - walkLength, softclipLength)),
                    reverseArray(Arrays.copyOfRange(refBases, refBases.length - walkLength, refBases.length)));
        }
        else
        {
            // a trailing softclip's M boundary is at its start, already boundary-outward
            int readStart = readBases.length - softclipLength;
            extended = maxScoringPrefix(
                    Arrays.copyOfRange(readBases, readStart, readStart + walkLength),
                    Arrays.copyOfRange(refBases, 0, walkLength));
        }

        if(extended == 0)
        {
            return 0;
        }

        applyExtension(working, softclipIndex, matchedIndex, extended);
        return extended;
    }

    private static void applyExtension(
            final List<CigarElement> working, final int softclipIndex, final int matchedIndex, final int extended)
    {
        CigarElement softclip = working.get(softclipIndex);
        CigarElement matched = working.get(matchedIndex);

        working.set(matchedIndex, new CigarElement(matched.getLength() + extended, matched.getOperator()));

        if(softclip.getLength() - extended == 0)
        {
            working.remove(softclipIndex);
        }
        else
        {
            working.set(softclipIndex, new CigarElement(softclip.getLength() - extended, CigarOperator.S));
        }
    }

    private static boolean isMatch(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ;
    }

    private static boolean isMatchedOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    private static int countIntrons(final List<CigarElement> elements)
    {
        int count = 0;
        for(CigarElement element : elements)
        {
            if(element.getOperator() == CigarOperator.N)
            {
                ++count;
            }
        }
        return count;
    }

    // one base pair's contribution to a bwa-mem score
    private static int baseScore(final byte readBase, final byte refBase)
    {
        return basesEqual(readBase, refBase) ? MATCH : MISMATCH;
    }

    // Longest prefix reaching the highest cumulative bwa-mem score (0 if none positive); the >= tie extends past an
    // internal mismatch.
    static int maxScoringPrefix(final byte[] read, final byte[] ref)
    {
        int length = Math.min(read.length, ref.length);
        int score = 0;
        int bestScore = 0;
        int bestLength = 0;
        for(int i = 0; i < length; ++i)
        {
            score += baseScore(read[i], ref[i]);
            if(score >= bestScore && score > 0)
            {
                bestScore = score;
                bestLength = i + 1;
            }
        }
        return bestLength;
    }

    // Straight bwa-mem score of read vs ref over their overlapping length.
    static int score(final byte[] read, final byte[] ref)
    {
        int total = 0;
        int length = Math.min(read.length, ref.length);
        for(int i = 0; i < length; ++i)
        {
            total += baseScore(read[i], ref[i]);
        }
        return total;
    }

    // bwa-mem-style genomic score of a lifted alignment: match/mismatch over aligned bases, affine gap per I/D, N and
    // S/H free.
    public int genomicScore(
            final String chromosome, final int alignmentStart, final String cigarStr, final byte[] readBases)
    {
        if(mRefGenome == null || readBases == null || cigarStr == null)
        {
            return Integer.MIN_VALUE;
        }

        Cigar cigar = CigarUtils.cigarFromStr(cigarStr);
        int totalScore = 0;
        int queryPos = 0;
        int refPos = alignmentStart;

        for(CigarElement element : cigar.getCigarElements())
        {
            int length = element.getLength();
            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    byte[] ref = refBases(chromosome, refPos, refPos + length - 1);
                    for(int i = 0; i < length; ++i)
                    {
                        // no reference, or a cigar reaching past SEQ (a mirrored failed-supp): score as a mismatch
                        boolean scorable = ref != null && queryPos + i < readBases.length;
                        totalScore += scorable ? baseScore(readBases[queryPos + i], ref[i]) : MISMATCH;
                    }
                    queryPos += length;
                    refPos += length;
                    break;

                case I:
                    totalScore += GAP_OPEN + (length - 1) * GAP_EXTEND;
                    queryPos += length;
                    break;

                case D:
                    totalScore += GAP_OPEN + (length - 1) * GAP_EXTEND;
                    refPos += length;
                    break;

                case N:
                    refPos += length;
                    break;

                case S:
                    queryPos += length;
                    break;

                default: // H, P consume neither query nor reference
                    break;
            }
        }

        return totalScore;
    }
}
