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

        List<SomaticCodingCount> count = SomaticCodingCount.create(winningAlleles);
        assertEquals(6, count.size());

        /* TODO - complete class impl first
        count = count.addVariant(true, CodingEffect.MISSENSE, Sets.newHashSet(HlaAllele.fromString("A*01:01"), HlaAllele.fromString("B*01:01")));
        count = count.addVariant(false, CodingEffect.MISSENSE, Sets.newHashSet(HlaAllele.fromString("B*01:01")));

        assertEquals(0.5, count.get(0).Total, 0.001)
        assertEquals(0.5, count.get(0).InframeIndel, 0.001)

        assertEquals(0.0, count.get(1).Total, 0.001)
        assertEquals(1.5, count.get(2).Total, 0.001)
        assertEquals(0.5, count.get(2).InframeIndel, 0.001)
        assertEquals(1.0, count.get(2).Missense, 0.001)

         */
    }

}
