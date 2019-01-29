package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Fusion;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.ImmutableFusion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneFusionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @NotNull
    private static Fusion fusionTestData() {
        return ImmutableFusion.builder()
                .reportable(true)
                .knownType("3P-Prom")
                .primarySource("")
                .clusterId("")
                .clusterCount("")
                .resolvedType("")
                .svIdUp("9104868")
                .chrUp("7")
                .posUp("106696332")
                .orientUp("1")
                .typeUp("DEL")
                .ploidyUp(1.65)
                .geneUp("PRKAR2B")
                .chrBandUp("q22.3")
                .transcriptUp("ENST00000265717")
                .strandUp("1")
                .regionTypeUp("Intronic")
                .codingTypeUp("Coding")
                .exonUp(1)
                .phaseUp("1")
                .exonMaxUp("11")
                .disruptiveUp("true")
                .exactBaseUp("-1")
                .codingBasesUp("307")
                .totalCodingUp("1254")
                .codingStartUp("106685353")
                .codingEndUp("106800027")
                .transStartUp("106685094")
                .transEndUp("106802256")
                .distancePrevUp("10673")
                .canonicalUp("true")
                .biotypeUp("protein_coding")
                .svIdDown("15444932")
                .chrDown("7")
                .posDown("116411792")
                .orientDown("-1")
                .typeDown("DEL")
                .ploidyDown(1.65)
                .geneDown("MET")
                .chrBandDown("q31.2")
                .transcriptDown("ENST00000318493")
                .strandDown("1")
                .regionTypeDown("Intronic")
                .codingTypeDown("Coding")
                .exonDown(14)
                .phaseDown("1")
                .exonMaxDown("21")
                .disruptiveDown("true")
                .exactBaseDown("-1")
                .codingBasesDown("1283")
                .totalCodingDown("4224")
                .codingStartDown("116339139")
                .codingEndDown("116436178")
                .transStartDown("116312459")
                .transEndDown("116436396")
                .distancePrevDown("84")
                .canonicalDown("true")
                .biotypeDown("protein_coding")
                .proteinsKept("Protein kinase domain")
                .proteinsLost("Sema domain")
                .build();
    }

    @Test
    public void canConvertGeneFusions() {
        List<ReportableGeneFusion> reportableFusions =
                ReportableGeneFusionFactory.fusionConvertToReportable(Lists.newArrayList(fusionTestData()));

        assertEquals(1, reportableFusions.size());

        ReportableGeneFusion reportableFusion = reportableFusions.get(0);
        assertEquals("PRKAR2B", reportableFusion.geneStart());
        assertEquals("MET", reportableFusion.geneEnd());
        assertEquals(1.65, reportableFusion.ploidy(), EPSILON);
    }
}