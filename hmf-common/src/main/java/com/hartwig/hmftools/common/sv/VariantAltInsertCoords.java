package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class VariantAltInsertCoords
{
    public final String Alt;
    public final String InsertSequence;
    public final Orientation Orient;

    public final String OtherChromsome;
    public final int OtherPosition;
    public final Orientation OtherOrient;

    public static final Pattern BREAKEND_REGEX = Pattern.compile("^(.*)([\\[\\]])(.+)[\\[\\]](.*)$");
    public static final Pattern SINGLE_BREAKEND_REGEX = Pattern.compile("^(([.].*)|(.*[.]))$");

    public static final char SINGLE_BREAKEND_CHAR = '.';
    public static final String SINGLE_BREAKEND_STR = String.valueOf(SINGLE_BREAKEND_CHAR);
    public static final byte SINGLE_BREAKEND_BYTE = (byte)SINGLE_BREAKEND_CHAR;

    public static final String POS_ORIENTATION_CHAR = "]";
    public static final String NEG_ORIENTATION_CHAR = "[";

    public VariantAltInsertCoords(
            final String alt, final String insertSequence, final Orientation orientation,
            final String otherChromsome, final int otherPosition, final Orientation otherOrientation)
    {
        Alt = alt;
        InsertSequence = insertSequence;
        Orient = orientation;
        OtherChromsome = otherChromsome;
        OtherPosition = otherPosition;
        OtherOrient = otherOrientation;
    }

    public static VariantAltInsertCoords fromRefAlt(final String altString, final String ref)
    {
        String alt = "";
        String insertSeq = "";
        Orientation orientation = null;

        String otherChromosome = "";
        int otherPosition = 0;
        Orientation otherOrientation = null;

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
                    orientation = FORWARD;
                }
                else
                {
                    String finalSequence = match.group(4);
                    insertSeq = !finalSequence.isEmpty() ? finalSequence.substring(0, finalSequence.length() - ref.length()) : "";
                    alt = altString.substring(altString.length() - 1);
                    orientation = REVERSE;
                }

                String orientStr = match.group(2);
                otherOrientation = orientStr.equals(POS_ORIENTATION_CHAR) ? FORWARD : REVERSE;

                String[] chrPos = match.group(3).split(":");
                otherChromosome = chrPos[0];
                otherPosition = Integer.parseInt(chrPos[1]);
            }
        }

        return new VariantAltInsertCoords(alt, insertSeq, orientation, otherChromosome, otherPosition, otherOrientation);
    }

    public static String formPairedAltString(
            final String alt, final String insertSequence, final String chromosome, int position, Orientation orientStart, Orientation orientEnd)
    {
        // pos orientation for DEL and DUP: AGAGATTATACTTTGTGTA[10:89712341[
        // pos orientation for INV: G]3:26664499]

        // neg orientation for DEL and DUP: ]10:89700299]GAGATTATACTTTGTGTAA
        // neg orientation for INV: [3:24566181[C

        if(orientStart.isForward() && orientEnd.isReverse())
            return String.format("%s%s[%s:%d[", alt, insertSequence, chromosome, position);
        else if(orientStart.isForward() && orientEnd.isForward())
            return String.format("%s%s]%s:%d]", alt, insertSequence, chromosome, position);
        else if(orientStart.isReverse() && orientEnd.isReverse())
            return String.format("[%s:%d[%s%s", chromosome, position, insertSequence, alt);
        else
            return String.format("]%s:%d]%s%s", chromosome, position, insertSequence, alt);
    }

    public static String formPairedAltString(
            final String alt, final String insertSequence, final String chromosome, int position, byte orientStart, byte orientEnd)
    {
        return formPairedAltString(
                alt, insertSequence, chromosome, position, Orientation.fromByte(orientStart), Orientation.fromByte(orientEnd));
    }

    public static String formSingleAltString(final String alt, final String insertSequence, Orientation orientation)
    {
        if(orientation.isForward())
            return String.format("%s%s.", alt, insertSequence);
        else
            return String.format(".%s%s", insertSequence, alt);
    }

    public static String formSingleAltString(final String alt, final String insertSequence, byte orientation)
    {
        return formSingleAltString(alt, insertSequence, Orientation.fromByte(orientation));
    }
}
