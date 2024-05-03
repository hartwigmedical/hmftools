package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.BREAKEND_REGEX;

import java.util.regex.Matcher;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class VariantAltInsertCoords
{
    public final String Alt;
    public final String InsertSequence;
    public final String Chromsome;
    public final int Position;
    public final Orientation Orient;

    public VariantAltInsertCoords(
            final String alt, final String insertSequence, final String chromsome, final int position, final Orientation orientation)
    {
        Alt = alt;
        InsertSequence = insertSequence;
        Chromsome = chromsome;
        Position = position;
        Orient = orientation;
    }

    public static VariantAltInsertCoords parseRefAlt(final String altString, final String ref)
    {
        String alt = "";
        String chromosome = "";
        String insertSeq = "";
        int position = 0;
        Orientation orientation = null;

        if(altString.startsWith("."))
        {
            alt = altString.substring(altString.length() - 1);
            insertSeq = altString.substring(ref.length(), altString.length() - 1);
            orientation = REVERSE;
        }
        else if(altString.endsWith("."))
        {
            alt = altString.substring(0, 1);
            insertSeq = altString.substring(1, altString.length() - ref.length());
            orientation = FORWARD;
        }
        else
        {
            final Matcher match = BREAKEND_REGEX.matcher(altString);

            if(match.matches())
            {
                if(match.group(1).length() > 0)
                {
                    String initialSequence = match.group(1);
                    insertSeq = !initialSequence.isEmpty() ? initialSequence.substring(ref.length()) : "";
                    alt = altString.substring(0, 1);
                }
                else
                {
                    String finalSequence = match.group(4);
                    insertSeq = !finalSequence.isEmpty() ? finalSequence.substring(0, finalSequence.length() - ref.length()) : "";
                    alt = altString.substring(altString.length() - 1);
                }

                String orientStr = match.group(2);
                orientation = orientStr.equals("]") ? FORWARD : REVERSE;

                String[] chrPos = match.group(3).split(":");
                chromosome = chrPos[0];
                position = Integer.parseInt(chrPos[1]);
            }
        }

        return new VariantAltInsertCoords(alt, insertSeq, chromosome, position, orientation);
    }

    public static String formPairedAltString(
            final String alt, final String insertSequence, final String chromosome, int position, Orientation orientStart, Orientation orientEnd)
    {
        if(orientStart.isForward() && orientEnd.isForward())
            return String.format("%s%s[%s:%d[", alt, insertSequence, chromosome, position);
        else if(orientStart.isForward() && orientEnd.isForward())
            return String.format("%s%s]%s:%d]", alt, insertSequence, chromosome, position);
        else if(orientStart.isReverse() && orientEnd.isReverse())
            return String.format("[%s:%d[%s%s", chromosome, position, insertSequence, alt);
        else
            return String.format("]%s:%d]%s%s", chromosome, position, insertSequence, alt);
    }

    public static String formSingleAltString(final String alt, final String insertSequence, Orientation orientation)
    {
        if(orientation.isForward())
            return String.format("%s%s.", alt, insertSequence);
        else
            return String.format(".%s%s", insertSequence, alt);
    }
}
