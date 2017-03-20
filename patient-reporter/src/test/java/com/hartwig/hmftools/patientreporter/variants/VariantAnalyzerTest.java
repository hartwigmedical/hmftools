package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.GenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class VariantAnalyzerTest {

    private static final String CHROMOSOME = "X";
    private static final String PASS_FILTER = "PASS";

    private static final String RIGHT_FEATURE_TYPE = "transcript";
    private static final String WRONG_FEATURE_TYPE = "sequence_feature";
    private static final String RIGHT_TRANSCRIPT = "TR";
    private static final String WRONG_TRANSCRIPT = "RT";
    private static final String REGION_ANNOTATION = RIGHT_TRANSCRIPT + ".1 (KODU)";

    @Test
    public void realCaseWorks() {
        final Slicer hmfSlicingRegion = SlicerFactory.fromSingleGenomeRegion(region(350, 450, REGION_ANNOTATION));
        final Slicer giabHighConfidenceRegion = SlicerFactory.fromSingleGenomeRegion(region(100, 1000));
        final Slicer cpctSlicingRegion = SlicerFactory.fromSingleGenomeRegion(region(400, 500));

        final VariantAnalyzer analyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion, giabHighConfidenceRegion,
                cpctSlicingRegion);

        final VariantAnnotation rightAnnotation = new VariantAnnotation.Builder().
                consequences(list(VariantConsequence.MISSENSE_VARIANT)).featureType(RIGHT_FEATURE_TYPE).
                featureID(RIGHT_TRANSCRIPT).build();

        final VariantAnnotation wrongTranscript = new VariantAnnotation.Builder().
                consequences(list(VariantConsequence.MISSENSE_VARIANT)).featureType(RIGHT_FEATURE_TYPE).
                featureID(WRONG_TRANSCRIPT).build();

        final VariantAnnotation wrongFeatureType = new VariantAnnotation.Builder().
                consequences(list(VariantConsequence.MISSENSE_VARIANT)).featureType(WRONG_FEATURE_TYPE).
                featureID(RIGHT_TRANSCRIPT).build();

        final VariantAnnotation wrongConsequence = new VariantAnnotation.Builder().
                consequences(list(VariantConsequence.OTHER)).featureType(RIGHT_FEATURE_TYPE).
                featureID(RIGHT_TRANSCRIPT).build();

        final List<SomaticVariant> variants = Lists.newArrayList(
                builder().position(420).annotations(Lists.newArrayList(rightAnnotation, wrongTranscript)).build(),
                builder().position(430).annotations(Lists.newArrayList(wrongConsequence)).build(),
                builder().position(440).annotations(Lists.newArrayList(wrongFeatureType)).build(),
                builder().position(460).annotations(Lists.newArrayList(rightAnnotation)).build());

        final VariantAnalysis analysis = analyzer.run(variants);

        assertEquals(4, analysis.allVariants().size());
        assertEquals(4, analysis.passedVariants().size());
        assertEquals(4, analysis.consensusPassedVariants().size());
        assertEquals(3, analysis.mutationalLoad());
        assertEquals(1, analysis.findings().size());
    }

    @NotNull
    private static List<VariantConsequence> list(@NotNull final VariantConsequence consequence) {
        return Lists.newArrayList(consequence);
    }

    @NotNull
    private static SomaticVariant.Builder builder() {
        return new SomaticVariant.Builder().type(VariantType.SNP).chromosome(CHROMOSOME).filter(PASS_FILTER);
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