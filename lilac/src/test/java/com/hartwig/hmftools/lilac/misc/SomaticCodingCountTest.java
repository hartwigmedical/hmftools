package com.hartwig.hmftools.lilac.misc;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

import org.junit.Test;

public class SomaticCodingCountTest
{
    @Test
    public void testRepeatAllele()
    {
        List<HlaAllele> winningAlleles = Lists.newArrayList(
                HlaAllele.fromString("A*01:01"), HlaAllele.fromString("A*01:01"),
                HlaAllele.fromString("B*01:01"), HlaAllele.fromString("B*01:02"),
                HlaAllele.fromString("C*01:01"), HlaAllele.fromString("C*01:03"));

        List<SomaticCodingCount> codingCounts = SomaticCodingCount.create(winningAlleles);
        assertEquals(6, codingCounts.size());

        SomaticCodingCount.addVariant(
                codingCounts, true, CodingEffect.MISSENSE,
                Lists.newArrayList(HlaAllele.fromString("A*01:01"), HlaAllele.fromString("B*01:01")));

        SomaticCodingCount.addVariant(codingCounts,
                false, CodingEffect.MISSENSE, Lists.newArrayList(HlaAllele.fromString("B*01:01")));

        assertEquals(0.5, codingCounts.get(0).total(), 0.001);
        assertEquals(0.5, codingCounts.get(0).inframeIndel(), 0.001);

        assertEquals(0.0, codingCounts.get(1).total(), 0.001);
        assertEquals(1.5, codingCounts.get(2).total(), 0.001);
        assertEquals(0.5, codingCounts.get(2).inframeIndel(), 0.001);
        assertEquals(1.0, codingCounts.get(2).missense(), 0.001);
    }

}
