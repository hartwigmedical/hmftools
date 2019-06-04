package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.structural.SvAnalysisDatamodelTestFactory.createTestFusionBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;

import org.junit.Test;

public class ReportableGeneFusionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canConvertGeneFusions() {
        List<ReportableGeneFusion> reportableFusions = Lists.newArrayList(
                createTestFusionBuilder().geneStart("PRKAR2B").geneEnd("MET").ploidy(1.65).build());

        assertEquals(1, reportableFusions.size());

        ReportableGeneFusion reportableFusion = reportableFusions.get(0);
        assertEquals("PRKAR2B", reportableFusion.geneStart());
        assertEquals("MET", reportableFusion.geneEnd());

        Double ploidy = reportableFusion.ploidy();
        assertNotNull(ploidy);
        assertEquals(1.65, ploidy, EPSILON);
    }
}