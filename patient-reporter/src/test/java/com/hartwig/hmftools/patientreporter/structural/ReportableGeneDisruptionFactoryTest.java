package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;
import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestDisruptionBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;

import org.junit.Test;

public class ReportableGeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertPairedDisruption() {
        ImmutableReportableDisruption.Builder pairedDisruptionBuilder =
                createTestDisruptionBuilder().svId(1).gene("ROPN1B").chromosome("3").chrBand("p12").type("INV").ploidy(1.12);

        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ROPN1B").minCopyNumber(1).maxCopyNumber(2).build());

        List<ReportableDisruption> pairedDisruptions = Lists.newArrayList(
                pairedDisruptionBuilder.exonUp(3).exonDown(4).build(),
                pairedDisruptionBuilder.exonUp(8).exonDown(9).build());

        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.convert(pairedDisruptions);

        assertEquals(1, reportableDisruptions.size());

        ReportableGeneDisruption disruption = reportableDisruptions.get(0);
        assertEquals("INV", disruption.type());
        assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals("Intron 3 -> Intron 8", disruption.range());
        assertEquals(3, disruption.firstAffectedExon());

        Double ploidy = disruption.ploidy();
        assertNotNull(ploidy);
        assertEquals(1.12, ploidy, EPSILON);
    }

    @Test
    public void doesNotPairDisruptionsOnDifferentGenes() {
        ImmutableReportableDisruption.Builder pairedDisruptionBuilder = createTestDisruptionBuilder().svId(1);

        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ROPN1B").minCopyNumber(1).maxCopyNumber(1).build(),
                        createTestCopyNumberBuilder().gene("SETD2").minCopyNumber(1).maxCopyNumber(1).build());

        List<ReportableDisruption> pairedDisruptions =
                Lists.newArrayList(pairedDisruptionBuilder.gene("ROPN1B").build(), pairedDisruptionBuilder.gene("SETD2").build());

        List<ReportableGeneDisruption> reportableDisruptions = ReportableGeneDisruptionFactory.convert(pairedDisruptions);

        assertEquals(2, reportableDisruptions.size());
    }
}