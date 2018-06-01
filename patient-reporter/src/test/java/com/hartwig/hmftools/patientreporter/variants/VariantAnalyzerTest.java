package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.snpeff.AnnotationTestFactory.createVariantAnnotationBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testMicrosatelliteAnalyzer;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantAnalyzerTest {

    private static final String PASS_FILTER = "PASS";

    private static final String RIGHT_FEATURE_TYPE = "transcript";
    private static final String WRONG_FEATURE_TYPE = "sequence_feature";
    private static final String RIGHT_TRANSCRIPT = "TR";
    private static final String WRONG_TRANSCRIPT = "RT";

    @Test
    public void realCaseWorks() {
        final VariantAnalyzer analyzer = VariantAnalyzer.of(Sets.newHashSet(RIGHT_TRANSCRIPT), testMicrosatelliteAnalyzer());

        final VariantAnnotation rightAnnotation =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(RIGHT_FEATURE_TYPE).
                        featureID(RIGHT_TRANSCRIPT).build();

        final VariantAnnotation wrongTranscript =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(RIGHT_FEATURE_TYPE)
                        .featureID(WRONG_TRANSCRIPT)
                        .build();

        final VariantAnnotation wrongFeatureType =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(WRONG_FEATURE_TYPE).
                        featureID(RIGHT_TRANSCRIPT).build();

        final VariantAnnotation wrongConsequence = createVariantAnnotationBuilder(VariantConsequence.OTHER).featureType(RIGHT_FEATURE_TYPE).
                featureID(RIGHT_TRANSCRIPT).build();

        final List<SomaticVariant> variants =
                Lists.newArrayList(builder().annotations(Lists.newArrayList(rightAnnotation, wrongTranscript)).build(),
                        builder().annotations(Lists.newArrayList(wrongTranscript)).build(),
                        builder().annotations(Lists.newArrayList(wrongConsequence)).build(),
                        builder().annotations(Lists.newArrayList(wrongFeatureType)).build(),
                        builder().annotations(Lists.newArrayList(rightAnnotation)).build());

        final VariantAnalysis analysis = analyzer.run(variants);

        assertEquals(5, analysis.passedVariants().size());
        assertEquals(0, analysis.mutationalLoad());
        assertEquals(2, analysis.variantReports().size());
    }

    @NotNull
    private static ImmutableSomaticVariantImpl.Builder builder() {
        return SomaticVariantTestBuilderFactory.create().filter(PASS_FILTER);
    }
}
