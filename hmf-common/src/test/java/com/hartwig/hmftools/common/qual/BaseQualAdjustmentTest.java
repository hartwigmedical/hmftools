package com.hartwig.hmftools.common.qual;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.phredQualToProbability;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.probabilityToPhredQualInt;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BaseQualAdjustmentTest
{
    @Test
    public void testBaseQualAdjustment()
    {
        double rawQual = 37;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        rawQual = 25;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        rawQual = 11;
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual + 1), 0.01);
        assertEquals(rawQual, BaseQualAdjustment.adjustBaseQual(rawQual - 1.49), 0.01);

        assertEquals(BASE_QUAL_MINIMUM, BaseQualAdjustment.adjustBaseQual(8.99));
    }

    @Test
    public void testQualCalcs()
    {
        byte qual = 25;
        double prob = phredQualToProbability(qual);
        assertEquals(0.003, prob, 0.001);

        double recalcQualDecimal = probabilityToPhredQual(prob);
        assertEquals(qual, recalcQualDecimal, 0.1);

        byte recalcQual = probabilityToPhredQualInt(prob);
        assertEquals(qual, recalcQual);
    }
}
