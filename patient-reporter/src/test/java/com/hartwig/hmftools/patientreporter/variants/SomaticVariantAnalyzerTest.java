package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.snpeff.AnnotationTestFactory.createVariantAnnotationBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testMicrosatelliteAnalyzer;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticVariantAnalyzerTest {

    private static final String PASS_FILTER = "PASS";

    private static final String RIGHT_FEATURE_TYPE = "transcript";
    private static final String WRONG_FEATURE_TYPE = "sequence_feature";
    private static final String RIGHT_TRANSCRIPT = "TR";
    private static final String WRONG_TRANSCRIPT = "RT";

    @Test
    public void realCaseWorks() {
        final SomaticVariantAnalyzer analyzer = SomaticVariantAnalyzer.of(Sets.newHashSet(RIGHT_TRANSCRIPT), testMicrosatelliteAnalyzer());

        final SnpEffAnnotation rightAnnotation =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(RIGHT_FEATURE_TYPE).
                        featureID(RIGHT_TRANSCRIPT).build();

        final SnpEffAnnotation wrongTranscript =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(RIGHT_FEATURE_TYPE)
                        .featureID(WRONG_TRANSCRIPT)
                        .build();

        final SnpEffAnnotation wrongFeatureType =
                createVariantAnnotationBuilder(VariantConsequence.MISSENSE_VARIANT).featureType(WRONG_FEATURE_TYPE).
                        featureID(RIGHT_TRANSCRIPT).build();

        final SnpEffAnnotation wrongConsequence = createVariantAnnotationBuilder(VariantConsequence.OTHER).featureType(RIGHT_FEATURE_TYPE).
                featureID(RIGHT_TRANSCRIPT).build();

        final List<EnrichedSomaticVariant> variants =
                Lists.newArrayList(builder().snpEffAnnotations(Lists.newArrayList(rightAnnotation, wrongTranscript)).build(),
                        builder().snpEffAnnotations(Lists.newArrayList(wrongTranscript)).build(),
                        builder().snpEffAnnotations(Lists.newArrayList(wrongConsequence)).build(),
                        builder().snpEffAnnotations(Lists.newArrayList(wrongFeatureType)).build(),
                        builder().snpEffAnnotations(Lists.newArrayList(rightAnnotation)).build());

        final SomaticVariantAnalysis analysis = analyzer.run(variants);

        assertEquals(5, analysis.passedVariants().size());
        assertEquals(0, analysis.mutationalLoad());
        assertEquals(2, analysis.variantReports().size());
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder builder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter(PASS_FILTER);
    }
}
