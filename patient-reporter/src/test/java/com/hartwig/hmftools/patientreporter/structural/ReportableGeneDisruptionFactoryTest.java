package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestCopyNumberBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Disruption;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.ImmutableDisruption;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneDisruptionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @NotNull
    private static Disruption disruptionTestData1() {
        return ImmutableDisruption.builder()
                .reportable(true)
                .svId("755779")
                .chromosome("3")
                .position("125593804")
                .orientation(-1)
                .type("INV")
                .ploidy(1)
                .gene("ROPN1B")
                .chrBand("p12")
                .transcript("ENST00000514116")
                .strand(1)
                .regionType("Upstream")
                .codingType("5P_UTR")
                .canonical(true)
                .biotype("protein_coding")
                .exonUp(0)
                .exonDown(1)
                .isDisruptive(false)
                .build();
    }

    @NotNull
    private static Disruption disruptionTestData2() {
        return ImmutableDisruption.builder()
                .reportable(true)
                .svId("9104967")
                .chromosome("3")
                .position("47079900")
                .orientation(1)
                .type("DUP")
                .ploidy(1)
                .gene("SETD2")
                .chrBand("p21.31")
                .transcript("ENST00000409792")
                .strand(-1)
                .regionType("Intronic")
                .codingType("Coding")
                .canonical(true)
                .biotype("protein_coding")
                .exonUp(17)
                .exonDown(18)
                .isDisruptive(true)
                .build();
    }

    @NotNull
    private static Disruption disruptionTestData3() {
        return ImmutableDisruption.builder()
                .reportable(true)
                .svId("9105167")
                .chromosome("22")
                .position("30041902")
                .orientation(-1)
                .type("DEL")
                .ploidy(1)
                .gene("NF2")
                .chrBand("q12.2")
                .transcript("ENST00000338641")
                .strand(1)
                .regionType("Intronic")
                .codingType("Coding")
                .canonical(true)
                .biotype("protein_coding")
                .exonUp(4)
                .exonDown(5)
                .isDisruptive(false)
                .build();
    }

    @Test
    public void canConvertPairedDisruption() {
        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ROPN1B").minCopyNumber(1).maxCopyNumber(1).build());
        List<Disruption> pairedDisruptions = Lists.newArrayList(disruptionTestData1());

        List<ReportableGeneDisruption> reportableDisruptions =
                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(pairedDisruptions, copyNumbers);

        assertEquals(1, reportableDisruptions.size());

        ReportableGeneDisruption disruption = reportableDisruptions.get(0);
        assertEquals("INV", disruption.type());
        assertEquals("3p12", disruption.location());
        assertEquals("ROPN1B", disruption.gene());
        assertEquals("Promoter Region Upstream", disruption.range());
        assertEquals(Integer.valueOf(1), disruption.geneMinCopies());
        assertEquals(Integer.valueOf(1), disruption.geneMaxCopies());
        assertEquals(0, disruption.firstAffectedExon());
        assertEquals(1D, disruption.ploidy(), EPSILON);
    }

    @Test
    public void doesNotPairDisruptionsOnDifferentGenes() {
        List<GeneCopyNumber> copyNumbers =
                Lists.newArrayList(createTestCopyNumberBuilder().gene("ROPN1B").minCopyNumber(1).maxCopyNumber(1).build(),
                        createTestCopyNumberBuilder().gene("SETD2").minCopyNumber(1).maxCopyNumber(1).build());
        List<Disruption> pairedDisruptions = Lists.newArrayList(disruptionTestData1(), disruptionTestData2());

        List<ReportableGeneDisruption> reportableDisruptions =
                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(pairedDisruptions, copyNumbers);

        assertEquals(2, reportableDisruptions.size());
    }

    @Test
    public void canConvertNormalDisruptionsWithoutCopyNumbers() {
        Disruption disruption1 = disruptionTestData1();
        Disruption disruption2 = disruptionTestData1();
        Disruption disruption3 = disruptionTestData1();

        List<Disruption> disruptions = Lists.newArrayList(disruption1, disruption2, disruption3);

        List<ReportableGeneDisruption> reportableDisruptions =
                ReportableGeneDisruptionFactory.disruptionConvertGeneDisruption(disruptions, Lists.newArrayList());
        assertEquals(0, reportableDisruptions.size());
    }
}