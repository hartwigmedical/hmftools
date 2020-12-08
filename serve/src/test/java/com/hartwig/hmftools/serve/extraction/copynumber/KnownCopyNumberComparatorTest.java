package com.hartwig.hmftools.serve.extraction.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class KnownCopyNumberComparatorTest {

    @Test
    public void canSortKnownCopyNumbers() {
        KnownCopyNumber copyNumber1 = ImmutableKnownCopyNumber.builder().gene("A").type(CopyNumberType.AMPLIFICATION).build();
        KnownCopyNumber copyNumber2 = ImmutableKnownCopyNumber.builder().gene("A").type(CopyNumberType.DELETION).build();
        KnownCopyNumber copyNumber3 = ImmutableKnownCopyNumber.builder().gene("B").type(CopyNumberType.AMPLIFICATION).build();
        KnownCopyNumber copyNumber4 = ImmutableKnownCopyNumber.builder().gene("C").type(CopyNumberType.DELETION).build();

        List<KnownCopyNumber> copyNumbers = Lists.newArrayList(copyNumber3, copyNumber2, copyNumber4, copyNumber1);
        copyNumbers.sort(new KnownCopyNumberComparator());

        assertEquals(copyNumber1, copyNumbers.get(0));
        assertEquals(copyNumber2, copyNumbers.get(1));
        assertEquals(copyNumber3, copyNumbers.get(2));
        assertEquals(copyNumber4, copyNumbers.get(3));
    }
}