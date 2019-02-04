package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestFusionBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class ReportableGeneFusionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertGeneFusions() {
        Fusion fusion = createTestFusionBuilder().geneUp("PRKAR2B").geneDown("MET").ploidyUp(1.65).ploidyDown(1.65).build();
        List<ReportableGeneFusion> reportableFusions =
                ReportableGeneFusionFactory.convert(Lists.newArrayList(fusion));

        assertEquals(1, reportableFusions.size());

        ReportableGeneFusion reportableFusion = reportableFusions.get(0);
        assertEquals("PRKAR2B", reportableFusion.geneStart());
        assertEquals("MET", reportableFusion.geneEnd());
        assertEquals(1.65, reportableFusion.ploidy(), EPSILON);
    }
}