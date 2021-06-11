package com.hartwig.hmftools.sage.rearrange;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RearrangementTest
{

    @Test
    public void testRearrangement1()
    {
        final String snv1 = ("          x      ");
        final String ref1 = ("CCT      TATGGTTT").replace(" ", "");
        final String alt1 = ("CCTGATGGGTCTGGTTT");
        final String ref2 = ("CCTTATGG      TTT").replace(" ", "");
        final String snv2 = ("   x             ");

        assertEquals(ref1, ref2);
        final List<VariantHotspot> victims = assertRearrangement(ref1, alt1, 7, 2, 7);

        // As inserts
        assertVariant(1002, "T", "TGATGGG", victims.get(0));
        assertVariant(1004, "A", "C", victims.get(1));
        assertVariant(1003, "T", "G", victims.get(2));
        assertVariant(1007, "G", "GGTCTGG", victims.get(3));

        // As deletes
        assertVariant(1002, "TGATGGG", "T", victims.get(4));
        assertVariant(1010, "C", "A", victims.get(5));
        assertVariant(1003, "G", "T", victims.get(6));
        assertVariant(1007, "GGTCTGG", "G", victims.get(7));
    }

    @Test
    public void testRearrangement2()
    {
        final String snv1 = ("          x      ");
        final String ref1 = ("CCT      TATGTTTT").replace(" ", "");
        final String alt1 = ("CCTGATGGGTCTGTTTT");
        final String ref2 = ("CCTTATG      TTTT").replace(" ", "");
        final String snv2 = ("   x             ");

        assertEquals(ref1, ref2);
        final List<VariantHotspot> victims = assertRearrangement(ref1, alt1, 7, 2, 6);

        // As inserts
        assertVariant(1002, "T", "TGATGGG", victims.get(0));
        assertVariant(1004, "A", "C", victims.get(1));
        assertVariant(1003, "T", "G", victims.get(2));
        assertVariant(1006, "G", "GGGTCTG", victims.get(3));

        // As deletes
        assertVariant(1002, "TGATGGG", "T", victims.get(4));
        assertVariant(1010, "C", "A", victims.get(5));
        assertVariant(1003, "G", "T", victims.get(6));
        assertVariant(1006, "GGGTCTG", "G", victims.get(7));
    }

    @Test
    public void testMnv1()
    {
        final String snv1 = ("          xx     ");
        final String ref1 = ("CCT      TATGGTTT").replace(" ", "");
        final String alt1 = ("CCTGATGGGTCAGGTTT");
        final String ref2 = ("CCTTATGG      TTT").replace(" ", "");
        final String snv2 = ("   x             ");

        assertEquals(ref1, ref2);
        final List<VariantHotspot> victims = assertRearrangement(ref1, alt1, 7, 2, 7);

        // As inserts
        assertVariant(1002, "T", "TGATGGG", victims.get(0));
        assertVariant(1004, "AT", "CA", victims.get(1));
        assertVariant(1003, "T", "G", victims.get(2));
        assertVariant(1007, "G", "GGTCAGG", victims.get(3));

        // As deletes
        assertVariant(1002, "TGATGGG", "T", victims.get(4));
        assertVariant(1010, "CA", "AT", victims.get(5));
        assertVariant(1003, "G", "T", victims.get(6));
        assertVariant(1007, "GGTCAGG", "G", victims.get(7));
    }

    @Test
    public void testMnv2()
    {
        final String snv1 = ("          x      ");
        final String ref1 = ("CCT      TATGGTTT").replace(" ", "");
        final String alt1 = ("CCTGTTGGGTCTGGTTT");
        final String ref2 = ("CCTTATGG      TTT").replace(" ", "");
        final String snv2 = ("   xx            ");

        assertEquals(ref1, ref2);
        final List<VariantHotspot> victims = assertRearrangement(ref1, alt1, 7, 2, 7);

        // As inserts
        assertVariant(1002, "T", "TGTTGGG", victims.get(0));
        assertVariant(1004, "A", "C", victims.get(1));
        assertVariant(1003, "TA", "GT", victims.get(2));
        assertVariant(1007, "G", "GGTCTGG", victims.get(3));

        // As deletes
        assertVariant(1002, "TGTTGGG", "T", victims.get(4));
        assertVariant(1010, "C", "A", victims.get(5));
        assertVariant(1003, "GT", "TA", victims.get(6));
        assertVariant(1007, "GGTCTGG", "G", victims.get(7));
    }

    @Test
    public void testLongIndelRearrangedMNVTooLong()
    {
        final String ref1 = ("CCT      TATGGTTT").replace(" ", "");
        final String alt1 = ("CCTCCCCCCTATGGTTT");

        final List<VariantHotspot> victims = Rearrangement.rearrangeInsert(1002, 7, 2, alt1.getBytes(), 2, ref1.getBytes());
        assertEquals(0, victims.size());
    }

    @Test
    public void testShortIndel()
    {
        final String ref1 = ("CCT  TATGGTTT").replace(" ", "");
        final String alt1 = ("CCTCCTCTGGTTT");

        final List<VariantHotspot> victims = Rearrangement.rearrangeInsert(1002, 3, 2, alt1.getBytes(), 2, ref1.getBytes());
        assertEquals(2, victims.size());
        assertVariant(1003, "TA", "CC", victims.get(0));
        assertVariant(1005, "T", "TCT", victims.get(1));
    }

    private List<VariantHotspot> assertRearrangement(String ref, String read, int indelLength, int indelIndexLeft, int indelIndexRight)
    {
        List<VariantHotspot> result = Lists.newArrayList();
        int refPosition = 1000;

        final List<VariantHotspot> leftInsert = Rearrangement.rearrangeInsert(refPosition + indelIndexLeft,
                indelLength,
                indelIndexLeft,
                read.getBytes(),
                indelIndexLeft,
                ref.getBytes());
        assertEquals(read, constructRead(refPosition, ref, leftInsert));
        assertEquals(refPosition + indelIndexRight, leftInsert.get(1).position());

        final List<VariantHotspot> leftDelete = Rearrangement.rearrangeDelete(refPosition + indelIndexLeft,
                indelLength,
                indelIndexLeft,
                ref.getBytes(),
                indelIndexLeft,
                read.getBytes());
        assertEquals(ref, constructRead(refPosition, read, leftDelete));
        assertEquals(refPosition + indelIndexRight, leftDelete.get(1).position());

        assertEquals(leftInsert.get(0).alt(), leftDelete.get(0).ref());
        assertEquals(leftInsert.get(0).ref(), leftDelete.get(0).alt());
        assertEquals(leftInsert.get(1).alt(), leftDelete.get(1).ref());
        assertEquals(leftInsert.get(1).ref(), leftDelete.get(1).alt());

        final List<VariantHotspot> rightInsert = Rearrangement.rearrangeInsert(refPosition + indelIndexRight,
                indelLength,
                indelIndexRight,
                read.getBytes(),
                indelIndexRight,
                ref.getBytes());
        assertEquals(read, constructRead(refPosition, ref, rightInsert));
        assertEquals(refPosition + indelIndexLeft, rightInsert.get(0).position());

        final List<VariantHotspot> rightDelete = Rearrangement.rearrangeDelete(refPosition + indelIndexRight,
                indelLength,
                indelIndexRight,
                ref.getBytes(),
                indelIndexRight,
                read.getBytes());
        assertEquals(ref, constructRead(refPosition, read, rightDelete));

        assertEquals(rightInsert.get(0).alt(), rightDelete.get(0).ref());
        assertEquals(rightInsert.get(0).ref(), rightDelete.get(0).alt());
        assertEquals(rightInsert.get(1).alt(), rightDelete.get(1).ref());
        assertEquals(rightInsert.get(1).ref(), rightDelete.get(1).alt());

        result.addAll(rightInsert);
        result.addAll(leftInsert);
        result.addAll(rightDelete);
        result.addAll(leftDelete);

        return result;
    }

    private String constructRead(int refPosition, @NotNull final String ref, @NotNull final List<VariantHotspot> variants)
    {
        int refIndex = 0;
        StringBuilder builder = new StringBuilder();

        for(VariantHotspot variant : variants)
        {
            int refLength = (int) variant.position() - refPosition;
            if(refLength > 0)
            {
                builder.append(ref, refIndex, refIndex + refLength);
                refIndex += refLength;
                refPosition += refLength;
            }

            builder.append(variant.alt());
            refIndex += variant.ref().length();
            refPosition += variant.ref().length();
        }

        builder.append(ref, refIndex, ref.length());

        return builder.toString();
    }

    private static void assertVariant(long position, String ref, String alt, VariantHotspot victim)
    {
        assertEquals(position, victim.position());
        assertEquals(ref, victim.ref());
        assertEquals(alt, victim.alt());
    }
}
