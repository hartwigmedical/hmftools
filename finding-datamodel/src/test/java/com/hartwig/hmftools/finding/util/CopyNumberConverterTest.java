package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class CopyNumberConverterTest  {

    private final static double EPSILON = 0.00001;

    @Test
    public void canFilterDisruptions() {
        Disruption disruption = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberConverter.filterDisruption(List.of(disruption)).get(0).reportedStatus());

        Disruption disruptionFilter = TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberConverter.filterDisruption(List.of(disruptionFilter)).get(0).reportedStatus());
    }

    @Test
    public void canFilterFusion() {
        Fusion fusion = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.1).build();
        assertEquals(ReportedStatus.REPORTED, CopyNumberConverter.filterFusion(List.of(fusion)).get(0).reportedStatus());

        Fusion fusionFilter = TestFindingFactory.fusionBuilder().junctionCopyNumber(0.03).build();
        assertEquals(ReportedStatus.CANDIDATE, CopyNumberConverter.filterFusion(List.of(fusionFilter)).get(0).reportedStatus());
    }

    @Test
    public void canRoundVariants() {
        List<SmallVariant> variant1 = CopyNumberConverter.convertSmallVariant(List.of(TestFindingFactory.variantBuilder().adjustedCopyNumber(0.144).build()));
        assertEquals(0.1, variant1.get(0).adjustedCopyNumber(), EPSILON);

        List<SmallVariant> variant2 = CopyNumberConverter.convertSmallVariant(List.of(TestFindingFactory.variantBuilder().adjustedCopyNumber(0.04).build()));
        assertEquals(0.0, variant2.get(0).adjustedCopyNumber(), EPSILON);
    }

    @Test
    public void canRoundDisruptions() {
        List<Disruption> disruption1 = CopyNumberConverter.convertDisruption(List.of(TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.144).build()));
        assertEquals(0.1, disruption1.get(0).disruptedCopyNumber(), EPSILON);

        List<Disruption> disruption2 = CopyNumberConverter.convertDisruption(List.of(TestFindingFactory.disruptionBuilder().disruptedCopyNumber(0.04).build()));
        assertEquals(0.0, disruption2.get(0).disruptedCopyNumber(), EPSILON);
    }

    @Test
    public void canRoundFusions() {
        List<Fusion> fusion1 = CopyNumberConverter.convertFusion(List.of(TestFindingFactory.fusionBuilder().junctionCopyNumber(0.144).build()));
        assertEquals(0.1, fusion1.get(0).junctionCopyNumber(), EPSILON);

        List<Fusion> fusion2 = CopyNumberConverter.convertFusion(List.of(TestFindingFactory.fusionBuilder().junctionCopyNumber(0.04).build()));
        assertEquals(0.0, fusion2.get(0).junctionCopyNumber(), EPSILON);

    }

    @Test
    public void canRoundHla() {
        List<HlaAllele> hla1 = CopyNumberConverter.convertHla(List.of(TestFindingFactory.hlaAlleleBuilder().tumorCopyNumber(0.144).build()));
        assertNotNull(hla1.get(0).tumorCopyNumber());
        assertEquals(0.1, hla1.get(0).tumorCopyNumber(), EPSILON);

        List<HlaAllele> hla2 = CopyNumberConverter.convertHla(List.of(TestFindingFactory.hlaAlleleBuilder().tumorCopyNumber(0.04).build()));
        assertNotNull(hla2.get(0).tumorCopyNumber());
        assertEquals(0.0, hla2.get(0).tumorCopyNumber(), EPSILON);
    }

    @Test
    public void testRoundCopyNumber() {
        assertEqualsDouble(2.0, CopyNumberConverter.roundCopyNumber(2.0136));
        assertEqualsDouble(2.5, CopyNumberConverter.roundCopyNumber(2.45));
        assertEqualsDouble(2.6, CopyNumberConverter.roundCopyNumber(2.61));
        assertEqualsDouble(0.0, CopyNumberConverter.roundCopyNumber(-2.61));
        assertEqualsDouble(114.0, CopyNumberConverter.roundCopyNumber(113.56));
        assertEqualsDouble(2.4, CopyNumberConverter.roundCopyNumber(2.44));
        assertNull(CopyNumberConverter.roundCopyNumberNullable(null));
    }

    public static void assertEqualsDouble(double expected, @Nullable Double actual) {
        assertNotNull(actual);
        assertEquals(expected, actual, EPSILON);
    }
}