package com.hartwig.hmftools.redux.splice.rescue;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// Small CIGAR helpers used by JunctionRescueResolver. Local to the rescue package so we can keep
// the rescue logic self-contained and testable without dragging in htsjdk's Cigar parser. The
// resolver only ever inspects boundary ops (leading/trailing S, the M chunks flanking softclip
// boundaries) so a regex parser is plenty.
public final class CigarShape
{
    public static final char OP_MATCH = 'M';
    public static final char OP_INSERTION = 'I';
    public static final char OP_DELETION = 'D';
    public static final char OP_SKIPPED = 'N';
    public static final char OP_SOFTCLIP = 'S';
    public static final char OP_HARDCLIP = 'H';
    public static final char OP_PAD = 'P';
    public static final char OP_SEQ_MATCH = '=';
    public static final char OP_SEQ_MISMATCH = 'X';

    private static final Pattern CIGAR_ELEMENT = Pattern.compile("(\\d+)([MIDNSHP=X])");

    public static class Element
    {
        public final int Length;
        public final char Op;

        public Element(final int length, final char op)
        {
            Length = length;
            Op = op;
        }

        @Override
        public String toString()
        {
            return Length + String.valueOf(Op);
        }
    }

    private CigarShape() {}

    public static List<Element> parse(final String cigar)
    {
        final List<Element> out = new ArrayList<>();
        if(cigar == null || cigar.isEmpty() || "*".equals(cigar))
            return out;

        final Matcher m = CIGAR_ELEMENT.matcher(cigar);
        int matchedTo = 0;
        while(m.find())
        {
            if(m.start() != matchedTo)
                throw new IllegalArgumentException("invalid CIGAR: " + cigar);
            out.add(new Element(Integer.parseInt(m.group(1)), m.group(2).charAt(0)));
            matchedTo = m.end();
        }
        if(matchedTo != cigar.length())
            throw new IllegalArgumentException("invalid CIGAR: " + cigar);
        return out;
    }

    public static String format(final List<Element> elements)
    {
        final StringBuilder sb = new StringBuilder();
        for(Element e : elements)
            sb.append(e.Length).append(e.Op);
        return sb.toString();
    }

    public static boolean consumesReference(final char op)
    {
        return op == OP_MATCH || op == OP_DELETION || op == OP_SKIPPED
                || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH;
    }

    public static boolean consumesRead(final char op)
    {
        return op == OP_MATCH || op == OP_INSERTION || op == OP_SOFTCLIP
                || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH;
    }

    // sum of read bases consumed by all ops (M/I/S/=/X). Used to validate that primary + supp
    // together cover the full read length with no overlap and no gap.
    public static int readLength(final List<Element> elements)
    {
        int total = 0;
        for(Element e : elements)
        {
            if(consumesRead(e.Op))
                total += e.Length;
        }
        return total;
    }

    // length of reference span (M + D + N + = + X). The position-1 + refSpan - 1 is the inclusive
    // end on the reference, used to derive intron start from primary and intron end from supp.
    public static int referenceSpan(final List<Element> elements)
    {
        int total = 0;
        for(Element e : elements)
        {
            if(consumesReference(e.Op))
                total += e.Length;
        }
        return total;
    }

    // read offset at which matched bases begin (skips a leading S if present, otherwise 0).
    public static int leadingSoftClip(final List<Element> elements)
    {
        if(elements.isEmpty() || elements.get(0).Op != OP_SOFTCLIP)
            return 0;
        return elements.get(0).Length;
    }

    public static int trailingSoftClip(final List<Element> elements)
    {
        if(elements.isEmpty())
            return 0;
        final Element last = elements.get(elements.size() - 1);
        return last.Op == OP_SOFTCLIP ? last.Length : 0;
    }

    public static boolean hasHardClip(final List<Element> elements)
    {
        for(Element e : elements)
        {
            if(e.Op == OP_HARDCLIP)
                return true;
        }
        return false;
    }

    // sum of M/=/X bases. The "matched" length under STAR's overhang sense (D/N don't count).
    public static int matchedBases(final List<Element> elements)
    {
        int total = 0;
        for(Element e : elements)
        {
            if(e.Op == OP_MATCH || e.Op == OP_SEQ_MATCH || e.Op == OP_SEQ_MISMATCH)
                total += e.Length;
        }
        return total;
    }
}
