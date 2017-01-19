package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.slicing.SlicerTestFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class VariantAnalyzerTest {

    private static final String CHROMOSOME = "X";

    @Test
    public void realCaseWorks() {
        final Slicer hmfSlicingRegion = SlicerTestFactory.forGenomeRegion(region(350, 450, "TRANS.1 (KODU)"));
        final Slicer giabHighConfidenceRegion = SlicerTestFactory.forGenomeRegion(region(100, 1000));
        final Slicer cpctSlicingRegion = SlicerTestFactory.forGenomeRegion(region(400, 500));

        final VariantAnalyzer analyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion,
                giabHighConfidenceRegion, cpctSlicingRegion);
        assertNotNull(analyzer);
    }

    @NotNull
    private static GenomeRegion region(final long start, final long end) {
        return region(start, end, null);
    }

    @NotNull
    private static GenomeRegion region(final long start, final long end, @Nullable String annotation) {
        return new GenomeRegion(CHROMOSOME, start, end, annotation);
    }
}