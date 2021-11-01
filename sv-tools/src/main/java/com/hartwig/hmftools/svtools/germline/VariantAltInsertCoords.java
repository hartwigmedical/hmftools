package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.BREAKEND_REGEX;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.regex.Matcher;

public class VariantAltInsertCoords
{
    public final String Alt;
    public final String InsertSequence;
    public final String Chromsome;
    public final int Position;
    public final byte Orientation;

    public VariantAltInsertCoords(
            final String alt, final String insertSequence, final String chromsome, final int position, final byte orientation)
    {
        Alt = alt;
        InsertSequence = insertSequence;
        Chromsome = chromsome;
        Position = position;
        Orientation = orientation;
    }

    public static VariantAltInsertCoords parseRefAlt(final String altString, final String ref)
    {
        String alt = "";
        String chromosome = "";
        String insertSeq = "";
        int position = 0;
        byte orientation = 0;

        if(altString.startsWith("."))
        {
            alt = altString.substring(altString.length() - 1);
            insertSeq = altString.substring(ref.length(), altString.length() - 1);
        }
        else if(altString.endsWith("."))
        {
            alt = altString.substring(0, 1);
            insertSeq = altString.substring(1, altString.length() - ref.length());
        }
        else
        {
            final Matcher match = BREAKEND_REGEX.matcher(altString);

            if(match.matches())
            {
                if(match.group(1).length() > 0)
                {
                    String initialSequence = match.group(1).substring(1);
                    insertSeq = initialSequence.substring(ref.length());
                    alt = altString.substring(0, 1);
                }
                else
                {
                    String finalSequence = match.group(4).substring(0, match.group(4).length() - 1);
                    insertSeq = finalSequence.substring(0, finalSequence.length() - ref.length());
                    alt = altString.substring(altString.length() - 1);
                }

                orientation = match.group(2).equals("]") ? POS_ORIENT : NEG_ORIENT;

                String[] chrPos = match.group(3).split(":");
                chromosome = chrPos[0];
                position = Integer.parseInt(chrPos[1]);
            }
        }

        return new VariantAltInsertCoords(alt, insertSeq, chromosome, position, orientation);
    }

    public static String formPairedAltString(
            final String alt, final String insertSequence, final String chromosome, int position, byte orientStart, byte orientEnd)
    {
        if(orientStart == POS_ORIENT && orientEnd == NEG_ORIENT)
            return String.format("%s%s[%s:%d[", alt, insertSequence, chromosome, position);
        else if(orientStart == POS_ORIENT && orientEnd == POS_ORIENT)
            return String.format("%s%s]%s:%d]", alt, insertSequence, chromosome, position);
        else if(orientStart == NEG_ORIENT && orientEnd == NEG_ORIENT)
            return String.format("[%s:%d[%s%s", chromosome, position, insertSequence, alt);
        else
            return String.format("]%s:%d]%s%s", chromosome, position, insertSequence, alt);
    }

    public static String formSingleAltString(final String alt, final String insertSequence, byte orientation)
    {
        if(orientation == POS_ORIENT)
            return String.format("%s%s.", alt, insertSequence);
        else
            return String.format(".%s%s", insertSequence, alt);
    }
}
