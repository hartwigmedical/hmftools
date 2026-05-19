package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class CopyNumberRoundingTransformerTest
{

    private final static double EPSILON = 0.00001;

    @Test
    public void canFilterDisruptions()
    {
        Disruption disruption = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberRoundingTransformer.filter(disruption).reportedStatus());

        Disruption disruptionFilter = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberRoundingTransformer.filter(disruptionFilter).reportedStatus());
    }

    @Test
    public void canFilterFusion()
    {
        Fusion fusion = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberRoundingTransformer.filter(fusion).reportedStatus());

        Fusion fusionFilter = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberRoundingTransformer.filter(fusionFilter).reportedStatus());
    }

    @Test
    public void canRoundVariants()
    {
        SmallVariant variant1 = CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.variantBuilder()
                .adjustedCopyNumber(0.144)
                .build());
        assertEquals(0.1, variant1.adjustedCopyNumber(), EPSILON);

        SmallVariant variant2 = CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.variantBuilder()
                .adjustedCopyNumber(0.04)
                .build());
        assertEquals(0.0, variant2.adjustedCopyNumber(), EPSILON);
    }

    @Test
    public void canRoundDisruptions()
    {
        Disruption disruption1 = CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.disruptionBuilder()
                .disruptedCopyNumber(0.144)
                .build());
        assertEquals(0.1, disruption1.disruptedCopyNumber(), EPSILON);

        Disruption disruption2 = CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.disruptionBuilder()
                .disruptedCopyNumber(0.04)
                .build());
        assertEquals(0.0, disruption2.disruptedCopyNumber(), EPSILON);
    }

    @Test
    public void canRoundFusions()
    {
        Fusion fusion1 = CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.fusionBuilder()
                .junctionCopyNumber(0.144)
                .build());
        assertEquals(0.1, fusion1.junctionCopyNumber(), EPSILON);

        Fusion fusion2 =
                CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.fusionBuilder().junctionCopyNumber(0.04).build());
        assertEquals(0.0, fusion2.junctionCopyNumber(), EPSILON);

    }

    @Test
    public void canRoundHla()
    {
        HlaAllele hla1 =
                CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.hlaAlleleBuilder().tumorCopyNumber(0.144).build());
        assertNotNull(hla1.tumorCopyNumber());
        assertEquals(0.1, hla1.tumorCopyNumber(), EPSILON);

        HlaAllele hla2 =
                CopyNumberRoundingTransformer.roundCopyNumbers(TestFindingFactory.hlaAlleleBuilder().tumorCopyNumber(0.04).build());
        assertNotNull(hla2.tumorCopyNumber());
        assertEquals(0.0, hla2.tumorCopyNumber(), EPSILON);
    }

    @Test
    public void testRoundCopyNumber()
    {
        assertEqualsDouble(2.0, CopyNumberRoundingTransformer.roundCopyNumber(2.0136));
        assertEqualsDouble(2.5, CopyNumberRoundingTransformer.roundCopyNumber(2.45));
        assertEqualsDouble(2.6, CopyNumberRoundingTransformer.roundCopyNumber(2.61));
        assertEqualsDouble(0.0, CopyNumberRoundingTransformer.roundCopyNumber(-2.61));
        assertEqualsDouble(114.0, CopyNumberRoundingTransformer.roundCopyNumber(113.56));
        assertEqualsDouble(2.4, CopyNumberRoundingTransformer.roundCopyNumber(2.44));
        assertNull(CopyNumberRoundingTransformer.roundCopyNumberNullable(null));
    }

    public static void assertEqualsDouble(double expected, @Nullable Double actual)
    {
        assertNotNull(actual);
        assertEquals(expected, actual, EPSILON);
    }
}