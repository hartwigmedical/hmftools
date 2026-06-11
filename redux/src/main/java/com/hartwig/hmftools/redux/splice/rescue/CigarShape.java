package com.hartwig.hmftools.redux.splice.rescue;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// CIGAR helpers for JunctionRescueResolver. Kept local to avoid htsjdk dependency in tests; a regex
// parser is sufficient since the resolver only inspects boundary ops.
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

    // Sum of read bases (M/I/S/=/X). Used to validate primary + supp cover the full read with no gap or overlap.
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

    // Reference span (M+D+N+=/X). pos-1 + refSpan - 1 gives the inclusive end; used to derive intron boundaries.
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

    // Leading soft-clip length, or 0 if none.
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

    // Sum of M/=/X bases (D/N excluded, matching STAR's overhang definition).
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
