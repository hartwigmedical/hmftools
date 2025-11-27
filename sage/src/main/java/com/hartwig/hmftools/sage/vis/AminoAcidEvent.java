package com.hartwig.hmftools.sage.vis;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER;

import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.tuple.Pair;

public interface AminoAcidEvent
{
    String FRAME_SHIFT = "fs";
    String TER = "Ter";
    char START_LOST = '?';
    char STOP = '*';
    String START_LOST_STR = "?";
    String STOP_STR = "*";

    int aminoAcidStartPos();
    int priority();

    default BaseRegion region() { return new BaseRegion(aminoAcidStartPos(), aminoAcidStartPos()); }

    record AminoAcidSubstitution(int aminoAcidStartPos, char alt) implements AminoAcidEvent
    {
        @Override
        public int priority()
        {
            return 4;
        }
    }

    record AminoAcidIns(int aminoAcidStartPos, int length) implements AminoAcidEvent
    {
        @Override
        public int priority()
        {
            return 5;
        }
    }

    record AminoAcidDel(int aminoAcidStartPos, int aminoAcidEndPos) implements AminoAcidEvent
    {
        @Override
        public int priority()
        {
            return 3;
        }

        @Override
        public BaseRegion region()
        {
            return new BaseRegion(aminoAcidStartPos, aminoAcidEndPos);
        }
    }

    record AminoAcidExt(int aminoAcidStartPos) implements AminoAcidEvent
    {
        @Override
        public int priority()
        {
            return 2;
        }
    }

    record AminoAcidStartLost(int aminoAcidStartPos) implements AminoAcidEvent
    {
        @Override
        public int priority()
        {
            return 1;
        }
    }

    static int compare(final AminoAcidEvent x, final AminoAcidEvent y)
    {
        if(x.aminoAcidStartPos() != y.aminoAcidStartPos())
            return x.aminoAcidStartPos() - y.aminoAcidStartPos();

        return x.priority() - y.priority();
    }

    static List<AminoAcidEvent> parse(final String hgvsStr)
    {
        CharSequence hgvs = new StringBuilder(hgvsStr);
        if(!"p.".contentEquals(hgvs.subSequence(0, 2)))
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        hgvs = hgvs.subSequence(2, hgvs.length());

        Pattern pattern = Pattern.compile("^([A-Z][a-z][a-z])(\\d+)(|[^\\d].*)$");
        Matcher matcher = pattern.matcher(hgvs);
        if(!matcher.find())
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        int matchLength = matcher.group(1).length() + matcher.group(2).length();
        String triAA = matcher.group(1);
        int aaPos = Integer.parseInt(matcher.group(2));
        Pair<String, Integer> ref1 = Pair.of(triAA, aaPos);
        hgvs = hgvs.subSequence(matchLength, hgvs.length());

        Pair<String, Integer> ref2 = null;
        if(!hgvs.isEmpty() && hgvs.charAt(0) == '_')
        {
            if(ref1.getLeft().equals(TER))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            hgvs = hgvs.subSequence(1, hgvs.length());
            matcher = pattern.matcher(hgvs);
            if(!matcher.find())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            matchLength = matcher.group(1).length() + matcher.group(2).length();
            triAA = matcher.group(1);
            aaPos = Integer.parseInt(matcher.group(2));
            ref2 = Pair.of(triAA, aaPos);
            hgvs = hgvs.subSequence(matchLength, hgvs.length());

            if(ref2.getLeft().equals(TER))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            if(ref2.getRight() <= ref1.getRight())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));
        }

        // synonymous
        if("=".contentEquals(hgvs))
            return Collections.emptyList();

        // frame shift
        if(FRAME_SHIFT.contentEquals(hgvs))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidExt(ref1.getRight()));
        }

        // stop lost
        if(TER.equals(ref1.getLeft()))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            pattern = Pattern.compile("^([A-Z][a-z][a-z])(ext\\*\\?)$");
            matcher = pattern.matcher(hgvs);
            if(!matcher.find())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            triAA = matcher.group(1);
            if(TER.equals(triAA))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidExt(ref1.getRight()));
        }

        if("ext*?".contentEquals(hgvs))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidExt(ref1.getRight()));
        }

        // simple missense and stop gained
        pattern = Pattern.compile("^([A-Z][a-z][a-z]|\\*)$");
        matcher = pattern.matcher(hgvs);
        if(matcher.find())
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            String triAltAA = hgvs.toString();
            char altAA = STOP_STR.equals(triAltAA) ? STOP : TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triAltAA).charAt(0);
            return List.of(new AminoAcidSubstitution(ref1.getRight(), altAA));
        }

        // start lost
        if(START_LOST_STR.contentEquals(hgvs))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidStartLost(ref1.getRight()));
        }

        // del
        if("del".contentEquals(hgvs))
        {
            int startPos = ref1.getRight();
            int endPos = ref2 == null ? startPos : ref2.getRight();
            return List.of(new AminoAcidDel(startPos, endPos));
        }

        // dup
        if("dup".contentEquals(hgvs))
        {
            int startPos = ref1.getRight();
            int endPos = ref2 == null ? startPos : ref2.getRight();
            int dupLength = endPos - startPos + 1;
            return List.of(new AminoAcidIns(endPos, dupLength));
        }

        if(!("ins".contentEquals(hgvs.subSequence(0, 3)) || "delins".contentEquals(hgvs.subSequence(0, 6))))
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        // ins, del-ins, including stop gained
        boolean isInsert = "ins".contentEquals(hgvs.subSequence(0, 3));
        if(isInsert)
            hgvs = hgvs.subSequence(3, hgvs.length());
        else
            hgvs = hgvs.subSequence(6, hgvs.length());

        if(ref2 == null)
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        if(isInsert && ref2.getRight() != ref1.getRight() + 1)
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        StringBuilder altAcids = new StringBuilder();
        pattern = Pattern.compile("^[A-Z][a-z][a-z].*$");
        while(!hgvs.isEmpty())
        {
            matcher = pattern.matcher(hgvs);
            if(matcher.find())
            {
                String triAltAA = hgvs.subSequence(0, 3).toString();
                if(TER.equals(triAltAA))
                    throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

                altAcids.append(TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triAltAA));
                hgvs = hgvs.subSequence(3, hgvs.length());
                continue;
            }

            if(STOP_STR.contentEquals(hgvs))
            {
                altAcids.append(STOP_STR);
                break;
            }

            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));
        }

        if(isInsert)
            return List.of(new AminoAcidIns(ref1.getRight(), altAcids.length()));

        List<AminoAcidEvent> events = Lists.newArrayList();
        int gapSize = ref2.getRight() - ref1.getRight() + 1;
        if(altAcids.length() > gapSize)
        {
            int insLen = altAcids.length() - gapSize;
            events.add(new AminoAcidIns(ref1.getRight() - 1, insLen));
            for(int i = ref1.getRight(); i <= ref2.getRight(); i++)
                events.add(new AminoAcidSubstitution(i, altAcids.charAt(i - ref1.getRight() + insLen)));

            return events;
        }

        int delLen = gapSize - altAcids.length();
        if(delLen > 0)
            events.add(new AminoAcidDel(ref1.getRight(), ref1.getRight() + delLen - 1));

        for(int i = 0; i < altAcids.length(); i++)
        {
            int pos = ref1.getRight() + i + delLen;
            events.add(new AminoAcidSubstitution(pos, altAcids.charAt(i)));
        }

        return events;
    }
}
