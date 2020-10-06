package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestDisruptionBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;

import org.junit.Test;

public class ReportableGeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertPairedDisruption() {
        ImmutableLinxBreakend.Builder pairedDisruptionBuilder =
                createTestDisruptionBuilder().svId(1).gene("ROPN1B").chromosome("3").chrBand("p12").type("INV").junctionCopyNumber(1.12);

        List<LinxBreakend> pairedDisruptions =
                Lists.newArrayList(pairedDisruptionBuilder.exonUp(3).exonDown(4).undisruptedCopyNumber(4.3).build(),
                        pairedDisruptionBuilder.exonUp(8).exonDown(9).undisruptedCopyNumber(2.1).build());

        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.convert(pairedDisruptions);

        assertEquals(1, reportableDisruptions.size());

        ReportableGeneDisruption disruption = reportableDisruptions.get(0);
        assertEquals("INV", disruption.type());
        assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals("Intron 3 -> Intron 8", disruption.range());
        assertEquals(3, disruption.firstAffectedExon());
        assertEquals(2.1, disruption.undisruptedCopyNumber(), EPSILON);

        Double copyNumber = disruption.junctionCopyNumber();
        assertNotNull(copyNumber);
        assertEquals(1.12, copyNumber, EPSILON);
    }

    @Test
    public void doesNotPairDisruptionsOnDifferentGenes() {
        ImmutableLinxBreakend.Builder pairedDisruptionBuilder = createTestDisruptionBuilder().svId(1);

        List<LinxBreakend> pairedDisruptions = Lists.newArrayList(pairedDisruptionBuilder.gene("ROPN1B")
                        .svId(1)
                        .junctionCopyNumber(1.0)
                        .undisruptedCopyNumber(1.0)
                        .build(),
                pairedDisruptionBuilder.gene("SETD2").svId(1).junctionCopyNumber(1.0).undisruptedCopyNumber(2.3).build(),
                pairedDisruptionBuilder.gene("SETD2").svId(1).junctionCopyNumber(1.0).undisruptedCopyNumber(1.7).build());

        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.convert(pairedDisruptions);

        assertEquals(2, reportableDisruptions.size());
    }
}