package com.hartwig.hmftools.esvee.assembly;

import java.util.Arrays;

import com.hartwig.hmftools.esvee.sequence.Sequence;

public class SupportCheckerTest
{
    private static SupportChecker checker()
    {
        return new SupportChecker();

        // TestUtils.config(
        //                Map.of("max_mismatched_count_for_strong_support", "1",
        //                        "max_mismatched_count_for_weak_support", "2"))
    }

    private static Sequence sequence(final String bases)
    {
        final byte[] baseQuality = new byte[bases.length()];
        Arrays.fill(baseQuality, (byte) 37);

        return Sequence.fromBytes(bases.getBytes(), baseQuality);
    }

    /* CHASHA FIXME
    @Test
    public void supportsWithSingleEdit()
    {
        final Sequence left = sequence("ATCGAAATGGGTC");
        final Sequence right = sequence("TCGAAACGGGTC");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);
    }

    @Test
    public void doesNotSupportsWithTwoEdits()
    {
        final Sequence left = sequence("ATCGAAATGGGTC");
        final Sequence right = sequence("TCGAAACTGGTC");

        assertTrue(checker().StrongSupport.supports(left, right)).isFalse();
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isNull();
    }

    @Test
    public void supportsWithSingleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGTC");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);

        assertTrue(checker().StrongSupport.supports(right, left)).isTrue();
        assertTrue(checker().StrongSupport.supportIndex(right, left)).isEqualTo(-1);
    }

    @Test
    public void doesNotSupportWithDoubleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGGTC");

        assertTrue(checker().StrongSupport.supports(left, right)).isFalse();
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isNull();
    }

    @Test
    public void nonContradictsWithDoubleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGGTC");

        assertTrue(checker().WeakSupport.supports(left, right)).isTrue();
        assertTrue(checker().WeakSupport.supportIndex(left, right)).isEqualTo(1);
        assertTrue(checker().WeakSupport.supportIndex(right, left)).isEqualTo(-1);
    }

    @Test
    public void supportsWithSingleContraction()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGTC");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);
    }

    @Test
    public void supportsWithSubstitutionNearRepeat()
    {
        final Sequence left = sequence("ATCGATCGAAAAAAAAAAATCGGCTA");
        final Sequence right = sequence("ATCGATCAAAAAAAAAAAATCGGCTA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithSubstitutionNearTwoRepeat()
    {
        final Sequence left = sequence( "ATCGATCGATATATATATATATATATATATTCGGCTA");
        final Sequence right = sequence("ATCGATCTATATATATATATATATATATATTCGGCTA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearStart()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAACAAAAAAAAAAAAAAAATCGATCGA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatAtStart()
    {
        //0123
        //GCTAAAAAAAAAAAAAAAAAAATCGATCGA
        //   AACAAAAAAAAAAAAAAAATCGATCGA
        final Sequence left = sequence( "GCTAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("AACAAAAAAAAAAAAAAAATCGATCGA");

        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(3);
    }

    @Test
    public void supportsWithInterruptedRepeatNearEnd()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAACAATCGATCGA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearEndAndChangedLength()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAACAATCGATCGA");

        assertTrue(checker().WeakSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearMiddle()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAACAAAAAAAATCGATCGA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithShortRepeatThatLooksInterruptedStrong()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATAAACGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAAAATAAACGATCGA");

        assertTrue(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithShortRepeatThatLooksInterruptedWeak()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATAAACGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAAATAAACGATCGA");

        assertTrue(checker().WeakSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsStartingInRepeat()
    {
        //32109876543210
        //             AAATAAAAAAAAA
        //GTCTGCAAAAAAAAAATA
        final Sequence left = sequence("AAATAAAAAAAAA");
        final Sequence right = sequence("GTCTGCAAAAAAAAAATA");

        assertTrue(checker().WeakSupport.supportIndex(left, right, 3)).isEqualTo(-13);
    }

    @Test
    public void supportsStartingInRepeatAmbiguous()
    {
        // There is "strong support" at index 0 and index 2, but obviously the one at index 2 is better.

        //012
        //ATATATATATATATATATGCTA
        //  ATATATATATATATATGCTA
        final Sequence left = sequence("ATATATATATATATATATGCTA");
        final Sequence right = sequence("ATATATATATATATATGCTA");

        assertTrue(checker().StrongSupport.supportsAt(left, right, 2)).isEqualTo(true);
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(2);
    }

    @Test
    public void supportsStartingInRepeatAmbiguousReverse()
    {
        // There is "strong support" at index 0 and index 2, but obviously the one at index 2 is better.

        //210
        //  ATATATATATATATATATGCTA
        //ATATATATATATATATATATGCTA
        final Sequence left = sequence("ATATATATATATATATATGCTA");
        final Sequence right = sequence("ATATATATATATATATATATGCTA");

        assertTrue(checker().StrongSupport.supportsAt(left, right, -2)).isEqualTo(true);
        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(-2);
    }

    @Test
    public void supportsStartingInRepeatPhaseChange()
    {
        //0123
        //TATATATATATATATATATGCTA
        //   ATATATATATATATATGCTA
        final Sequence left = sequence("TATATATATATATATATATGCTA");
        final Sequence right = sequence("ATATATATATATATATGCTA");

        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(3);
    }

    @Test
    public void supportsStartingInRepeatReversePhaseChange()
    {
        //10
        // ATATATATATATATATATGCTA
        //TATATATATATATATATATGCTA
        final Sequence left = sequence("ATATATATATATATATATGCTA");
        final Sequence right = sequence("TATATATATATATATATATGCTA");

        assertTrue(checker().StrongSupport.supportIndex(left, right)).isEqualTo(-1);
    }

    @Test
    public void supportsStartingInRepeatReversePhaseChange2()
    {
        //6543210
        //      ATATATATATATATATGCTATGCTTTGGCTGCTTTAC
        //GCTTATATATATATATATATATGCTATGCTTTGGCT
        final Sequence left = sequence("ATATATATATATATATGCTATGCTTTGGCTGCTTTAC");
        final Sequence right = sequence("GCTTATATATATATATATATATGCTATGCTTTGGCT");

        assertTrue(checker().StrongSupport.supportIndex(left, right, 20)).isEqualTo(-6);
    }

     */
}
