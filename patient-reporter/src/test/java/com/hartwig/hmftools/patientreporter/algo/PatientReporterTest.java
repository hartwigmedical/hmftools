package com.hartwig.hmftools.patientreporter.algo;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReporterData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testHmfReporterData;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.BaseReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.VariantAnnotator;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String FUSIONS_CSV = Resources.getResource("csv").getPath() + File.separator + "cosmic_gene_fusions.csv";

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfPanelGeneMap());
        final BaseReporterData baseReporterData = testBaseReporterData();
        final HmfReporterData reporterData = testHmfReporterData();
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.of(geneModel, reporterData.microsatelliteAnalyzer());
        final StructuralVariantAnalyzer svAnalyzer =
                new StructuralVariantAnalyzer(new TestAnnotator(), geneModel.regions(), CosmicFusions.readFromCSV(FUSIONS_CSV));
        final PatientReporter algo = ImmutablePatientReporter.of(baseReporterData, reporterData, variantAnalyzer, svAnalyzer);
        assertNotNull(algo.run(RUN_DIRECTORY, null));
    }

    private static class TestAnnotator implements VariantAnnotator {

        @Override
        @NotNull
        public List<StructuralVariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
            final List<StructuralVariantAnnotation> result = Lists.newArrayList();
            for (final StructuralVariant sv : variants) {
                final StructuralVariantAnnotation ann = new StructuralVariantAnnotation(sv);
                final GeneAnnotation g1 =
                        new GeneAnnotation(ann, true, "PNPLA7", Collections.singletonList("PNPLA7"), "ENSG00000130653", -1);
                g1.addTranscript(new Transcript(g1, "ENST00000406427", 12, 0, 13, 0, 37, true));
                ann.annotations().add(g1);

                final GeneAnnotation g2 =
                        new GeneAnnotation(ann, false, "TMPRSS2", Collections.singletonList("TMPRSS2"), "ENSG00000184012", -1);
                g2.addTranscript(new Transcript(g2, "ENST00000398585", 1, 0, 2, 0, 14, true));
                ann.annotations().add(g2);

                result.add(ann);
            }
            return result;
        }
    }
}