package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusions;
import com.hartwig.hmftools.patientreporter.variants.StructuralVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.Breakend;
import com.hartwig.hmftools.svannotation.GeneAnnotation;
import com.hartwig.hmftools.svannotation.Transcript;
import com.hartwig.hmftools.svannotation.VariantAnnotation;
import com.hartwig.hmftools.svannotation.VariantAnnotator;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String FUSIONS_CSV = Resources.getResource("csv").getPath() + File.separator + "cosmic_gene_fusions.csv";

    private static class TestAnnotator implements VariantAnnotator {

        @Override
        public List<VariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
            final List<VariantAnnotation> result = Lists.newArrayList();
            for (final StructuralVariant sv : variants) {
                final VariantAnnotation ann = new VariantAnnotation(sv);

                final Breakend b1 = new Breakend(ann, sv.startChromosome(), sv.startPosition(), sv.startOrientation(), sv.startAF());
                final GeneAnnotation g1 = new GeneAnnotation(b1, "PNPLA7", Collections.singletonList("PNPLA7"), "ENSG00000130653", -1);
                g1.addTranscriptAnnotation(new Transcript(g1, "ENST00000406427", 12, 0, 13, 0, 37, true));
                b1.addGeneAnnotation(g1);

                final Breakend b2 = new Breakend(ann, sv.endChromosome(), sv.endPosition(), sv.endOrientation(), sv.endAF());
                final GeneAnnotation g2 = new GeneAnnotation(b2, "TMPRSS2", Collections.singletonList("TMPRSS2"), "ENSG00000184012", -1);
                g2.addTranscriptAnnotation(new Transcript(g2, "ENST00000398585", 1, 0, 2, 0, 14, true));
                b2.addGeneAnnotation(g2);

                ann.setBreakendAnnotations(b1, b2);

                result.add(ann);
            }
            return result;
        }

        @Override
        public VariantAnnotation annotateRegion(final GenomeRegion region) {
            return null;
        }
    }

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException, DRException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.asMap());
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(geneModel);
        final StructuralVariantAnalyzer svAnalyzer =
                new StructuralVariantAnalyzer(new TestAnnotator(), geneModel.hmfRegions(), COSMICGeneFusions.readFromCSV(FUSIONS_CSV));
        final PatientReporter algo =
                new PatientReporter(buildTestCpctEcrfModel(), LimsJsonModel.buildEmptyModel(), variantAnalyzer, svAnalyzer);
        assertNotNull(algo.run(RUN_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildTestCpctEcrfModel() {
        return new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
    }
}