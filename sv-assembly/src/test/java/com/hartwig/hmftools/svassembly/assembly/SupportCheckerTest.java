package com.hartwig.hmftools.svassembly.assembly;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.Arrays;
import java.util.Map;

import com.hartwig.hmftools.svassembly.TestUtils;
import com.hartwig.hmftools.svassembly.models.Sequence;

import org.junit.Test;

public class SupportCheckerTest
{
    private static SupportChecker checker()
    {
        return new SupportChecker(TestUtils.config(
                Map.of("max_mismatched_count_for_strong_support", "1",
                        "max_mismatched_count_for_weak_support", "2")));
    }

    private static Sequence sequence(final String bases)
    {
        final byte[] baseQuality = new byte[bases.length()];
        Arrays.fill(baseQuality, (byte) 37);

        return Sequence.fromBytes(bases.getBytes(), baseQuality);
    }

    @Test
    public void supportsWithSingleEdit()
    {
        final Sequence left = sequence("ATCGAAATGGGTC");
        final Sequence right = sequence("TCGAAACGGGTC");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);
    }

    @Test
    public void doesNotSupportsWithTwoEdits()
    {
        final Sequence left = sequence("ATCGAAATGGGTC");
        final Sequence right = sequence("TCGAAACTGGTC");

        assertThat(checker().StrongSupport.supports(left, right)).isFalse();
        assertThat(checker().StrongSupport.supportIndex(left, right)).isNull();
    }

    @Test
    public void supportsWithSingleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGTC");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);

        assertThat(checker().StrongSupport.supports(right, left)).isTrue();
        assertThat(checker().StrongSupport.supportIndex(right, left)).isEqualTo(-1);
    }

    @Test
    public void doesNotSupportWithDoubleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGGTC");

        assertThat(checker().StrongSupport.supports(left, right)).isFalse();
        assertThat(checker().StrongSupport.supportIndex(left, right)).isNull();
    }

    @Test
    public void nonContradictsWithDoubleExtension()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGGGGTC");

        assertThat(checker().WeakSupport.supports(left, right)).isTrue();
        assertThat(checker().WeakSupport.supportIndex(left, right)).isEqualTo(1);
        assertThat(checker().WeakSupport.supportIndex(right, left)).isEqualTo(-1);
    }

    @Test
    public void supportsWithSingleContraction()
    {
        final Sequence left = sequence("ATCGAAATGGGGGGGGTC");
        final Sequence right = sequence("TCGAAATGGGGGGGTC");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(1);
    }

    @Test
    public void supportsWithSubstitutionNearRepeat()
    {
        final Sequence left = sequence("ATCGATCGAAAAAAAAAAATCGGCTA");
        final Sequence right = sequence("ATCGATCAAAAAAAAAAAATCGGCTA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithSubstitutionNearTwoRepeat()
    {
        final Sequence left = sequence( "ATCGATCGATATATATATATATATATATATTCGGCTA");
        final Sequence right = sequence("ATCGATCTATATATATATATATATATATATTCGGCTA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearStart()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAACAAAAAAAAAAAAAAAATCGATCGA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatAtStart()
    {
        //0123
        //GCTAAAAAAAAAAAAAAAAAAATCGATCGA
        //   AACAAAAAAAAAAAAAAAATCGATCGA
        final Sequence left = sequence( "GCTAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("AACAAAAAAAAAAAAAAAATCGATCGA");

        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(3);
    }

    @Test
    public void supportsWithInterruptedRepeatNearEnd()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAACAATCGATCGA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearEndAndChangedLength()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAACAATCGATCGA");

        assertThat(checker().WeakSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithInterruptedRepeatNearMiddle()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATCGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAACAAAAAAAATCGATCGA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithShortRepeatThatLooksInterruptedStrong()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATAAACGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAAAATAAACGATCGA");

        assertThat(checker().StrongSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsWithShortRepeatThatLooksInterruptedWeak()
    {
        final Sequence left = sequence( "ATCGGCTATCAAAAAAAAAAAAAAAAAAATAAACGATCGA");
        final Sequence right = sequence("ATCGGCTATCAAAAAAAAAAAAAAAAATAAACGATCGA");

        assertThat(checker().WeakSupport.supports(left, right)).isTrue();
    }

    @Test
    public void supportsStartingInRepeat()
    {
        //32109876543210
        //             AAATAAAAAAAAA
        //GTCTGCAAAAAAAAAATA
        final Sequence left = sequence("AAATAAAAAAAAA");
        final Sequence right = sequence("GTCTGCAAAAAAAAAATA");

        assertThat(checker().WeakSupport.supportIndex(left, right, 3)).isEqualTo(-13);
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

        assertThat(checker().StrongSupport.supportsAt(left, right, 2)).isEqualTo(true);
        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(2);
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

        assertThat(checker().StrongSupport.supportsAt(left, right, -2)).isEqualTo(true);
        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(-2);
    }

    @Test
    public void supportsStartingInRepeatPhaseChange()
    {
        //0123
        //TATATATATATATATATATGCTA
        //   ATATATATATATATATGCTA
        final Sequence left = sequence("TATATATATATATATATATGCTA");
        final Sequence right = sequence("ATATATATATATATATGCTA");

        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(3);
    }

    @Test
    public void supportsStartingInRepeatReversePhaseChange()
    {
        //10
        // ATATATATATATATATATGCTA
        //TATATATATATATATATATGCTA
        final Sequence left = sequence("ATATATATATATATATATGCTA");
        final Sequence right = sequence("TATATATATATATATATATGCTA");

        assertThat(checker().StrongSupport.supportIndex(left, right)).isEqualTo(-1);
    }

    @Test
    public void supportsStartingInRepeatReversePhaseChange2()
    {
        //6543210
        //      ATATATATATATATATGCTATGCTTTGGCTGCTTTAC
        //GCTTATATATATATATATATATGCTATGCTTTGGCT
        final Sequence left = sequence("ATATATATATATATATGCTATGCTTTGGCTGCTTTAC");
        final Sequence right = sequence("GCTTATATATATATATATATATGCTATGCTTTGGCT");

        assertThat(checker().StrongSupport.supportIndex(left, right, 20)).isEqualTo(-6);
    }
}
