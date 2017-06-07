package com.hartwig.hmftools.patientreporter.report;

import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotationFactory;

import org.junit.Test;

import net.sf.jasperreports.engine.JRDataSource;

public class GenePanelDataSourceTest {

    @Test
    public void canCreateGenePanelFor3Genes() {
        final List<HMFSlicingAnnotation> genes = Lists.newArrayList(HMFSlicingAnnotationFactory.create("XX", 1, "A"),
                HMFSlicingAnnotationFactory.create("XX", 2, "B"), HMFSlicingAnnotationFactory.create("XX", 3, "C"));

        final JRDataSource dataSource = GenePanelDataSource.fromHMFSlicingAnnotations(genes);
        assertNotNull(dataSource);
    }
}