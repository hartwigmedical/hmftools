package com.hartwig.hmftools.sage.vis;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_FRAMESHIFT;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_PROTEIN_ID;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_START_LOST;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_STOP_GAINED;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_STOP_LOST;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_TYPE_DEL;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_TYPE_DEL_INS;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_TYPE_DUP;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_TYPE_INS;

import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.tuple.Pair;

public interface AminoAcidEvent
{
    String TER = "Ter";
    char START_LOST_CHAR = HGVS_START_LOST.charAt(0);
    char STOP_CHAR = HGVS_STOP_GAINED.charAt(0);

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
        CharSequence hgvsWorking = new StringBuilder(hgvsStr);

        if(!HGVS_PROTEIN_ID.contentEquals(hgvsWorking.subSequence(0, 2)))
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        hgvsWorking = hgvsWorking.subSequence(2, hgvsWorking.length());

        Pattern pattern = Pattern.compile("^([A-Z][a-z][a-z])(\\d+)(|[^\\d].*)$");
        Matcher matcher = pattern.matcher(hgvsWorking);
        if(!matcher.find())
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        int matchLength = matcher.group(1).length() + matcher.group(2).length();
        String triAA = matcher.group(1);
        int aaPos = Integer.parseInt(matcher.group(2));

        // make a pair of amino acid character and index in the reference
        Pair<String, Integer> ref1 = Pair.of(triAA, aaPos);
        hgvsWorking = hgvsWorking.subSequence(matchLength, hgvsWorking.length());

        // make a second pair if the reference spans more than one amino acid
        Pair<String, Integer> ref2 = null;
        if(!hgvsWorking.isEmpty() && hgvsWorking.charAt(0) == '_')
        {
            if(ref1.getLeft().equals(TER))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            hgvsWorking = hgvsWorking.subSequence(1, hgvsWorking.length());
            matcher = pattern.matcher(hgvsWorking);
            if(!matcher.find())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            matchLength = matcher.group(1).length() + matcher.group(2).length();
            triAA = matcher.group(1);
            aaPos = Integer.parseInt(matcher.group(2));
            ref2 = Pair.of(triAA, aaPos);
            hgvsWorking = hgvsWorking.subSequence(matchLength, hgvsWorking.length());

            if(ref2.getLeft().equals(TER))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            if(ref2.getRight() <= ref1.getRight())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));
        }

        // synonymous
        if("=".contentEquals(hgvsWorking))
            return Collections.emptyList();

        // frame shift
        if(HGVS_FRAMESHIFT.contentEquals(hgvsWorking))
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
            matcher = pattern.matcher(hgvsWorking);
            if(!matcher.find())
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            triAA = matcher.group(1);
            if(TER.equals(triAA))
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidExt(ref1.getRight()));
        }

        if(HGVS_STOP_LOST.contentEquals(hgvsWorking))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidExt(ref1.getRight()));
        }

        // simple missense and stop gained
        pattern = Pattern.compile("^([A-Z][a-z][a-z]|\\*)$");
        matcher = pattern.matcher(hgvsWorking);
        if(matcher.find())
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            String triAltAA = hgvsWorking.toString();
            char altAA = HGVS_STOP_GAINED.equals(triAltAA) ? STOP_CHAR : TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triAltAA).charAt(0);
            return List.of(new AminoAcidSubstitution(ref1.getRight(), altAA));
        }

        // start lost
        if(HGVS_START_LOST.contentEquals(hgvsWorking))
        {
            if(ref2 != null)
                throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

            return List.of(new AminoAcidStartLost(ref1.getRight()));
        }

        // del
        if(HGVS_TYPE_DEL.contentEquals(hgvsWorking))
        {
            int startPos = ref1.getRight();
            int endPos = ref2 == null ? startPos : ref2.getRight();
            return List.of(new AminoAcidDel(startPos, endPos));
        }

        // dup
        if(HGVS_TYPE_DUP.contentEquals(hgvsWorking))
        {
            int startPos = ref1.getRight();
            int endPos = ref2 == null ? startPos : ref2.getRight();
            int dupLength = endPos - startPos + 1;
            return List.of(new AminoAcidIns(endPos, dupLength));
        }

        if(!(HGVS_TYPE_INS.contentEquals(hgvsWorking.subSequence(0, 3)) || HGVS_TYPE_DEL_INS.contentEquals(hgvsWorking.subSequence(0, 6))))
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        // ins, del-ins, including stop gained
        boolean isInsert = HGVS_TYPE_INS.contentEquals(hgvsWorking.subSequence(0, 3));
        boolean isDelInsert = !isInsert;

        if(isInsert)
            hgvsWorking = hgvsWorking.subSequence(3, hgvsWorking.length());
        else
            hgvsWorking = hgvsWorking.subSequence(6, hgvsWorking.length());

        if(ref2 == null && isDelInsert)
        {
            ref2 = Pair.of(hgvsWorking.toString(), ref1.getRight());
            // throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));
        }

        if(isInsert && ref2.getRight() != ref1.getRight() + 1)
            throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

        StringBuilder altAcids = new StringBuilder();
        pattern = Pattern.compile("^[A-Z][a-z][a-z].*$");
        while(!hgvsWorking.isEmpty())
        {
            matcher = pattern.matcher(hgvsWorking);
            if(matcher.find())
            {
                String triAltAA = hgvsWorking.subSequence(0, 3).toString();
                if(TER.equals(triAltAA))
                    throw new IllegalArgumentException(format("hgvs(%s) does not follow expected format", hgvsStr));

                altAcids.append(TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triAltAA));
                hgvsWorking = hgvsWorking.subSequence(3, hgvsWorking.length());
                continue;
            }

            if(HGVS_STOP_GAINED.contentEquals(hgvsWorking))
            {
                altAcids.append(HGVS_STOP_GAINED);
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
