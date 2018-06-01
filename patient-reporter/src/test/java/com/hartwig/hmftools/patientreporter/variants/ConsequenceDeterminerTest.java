package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.snpeff.AnnotationTestFactory.createVariantAnnotationBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.ImmutableVariantAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

import org.junit.Test;

public class ConsequenceDeterminerTest {

    private static final String CHROMOSOME = "X";
    private static final long POSITION = 150L;
    private static final String GENE = "GENE";

    private static final String TRANSCRIPT = "TRANS.1";

    private static final String REF = "R";
    private static final String ALT = "A";
    private static final String COSMIC_ID = "123";
    private static final int ALLELE_READ_COUNT = 1;
    private static final int TOTAL_READ_COUNT = 2;

    private static final String HGVS_CODING = "c.RtoA";
    private static final String HGVS_PROTEIN = "p.RtoA";

    @Test
    public void worksAsExpected() {
        final ConsequenceDeterminer determiner = new ConsequenceDeterminer(Sets.newHashSet(TRANSCRIPT));

        final VariantConsequence rightConsequence = VariantConsequence.MISSENSE_VARIANT;
        final VariantConsequence wrongConsequence = VariantConsequence.OTHER;

        final ImmutableVariantAnnotation.Builder annotationBuilder = createVariantAnnotationBuilder().featureID(TRANSCRIPT).
                featureType(ConsequenceDeterminer.FEATURE_TYPE_TRANSCRIPT).gene(GENE).hgvsCoding(HGVS_CODING).
                hgvsProtein(HGVS_PROTEIN);
        final VariantAnnotation rightAnnotation = annotationBuilder.consequences(Lists.newArrayList(rightConsequence)).build();
        final VariantAnnotation wrongAnnotation = annotationBuilder.consequences(Lists.newArrayList(wrongConsequence)).build();

        final ImmutableSomaticVariantImpl.Builder variantBuilder = SomaticVariantTestBuilderFactory.create().
                chromosome(CHROMOSOME).ref(REF).alt(ALT).cosmicID(COSMIC_ID).position(POSITION).
                totalReadCount(TOTAL_READ_COUNT).alleleReadCount(ALLELE_READ_COUNT);

        final SomaticVariant rightVariant = variantBuilder.annotations(Lists.newArrayList(rightAnnotation)).build();
        final SomaticVariant wrongVariant = variantBuilder.annotations(Lists.newArrayList(wrongAnnotation)).build();

        final List<VariantReport> variantReports = determiner.run(Lists.newArrayList(rightVariant, wrongVariant));
        assertEquals(1, variantReports.size());

        final VariantReport variantReport = variantReports.get(0);
        assertEquals(GENE, variantReport.gene());
        assertEquals(CHROMOSOME + ":" + POSITION, variantReport.variant().chromosomePosition());
        assertEquals(REF, variantReport.variant().ref());
        assertEquals(ALT, variantReport.variant().alt());
        assertEquals(TRANSCRIPT, variantReport.transcript());
        assertEquals(HGVS_CODING, variantReport.hgvsCoding());
        assertEquals(HGVS_PROTEIN, variantReport.hgvsProtein());
        assertEquals(rightConsequence.readableSequenceOntologyTerm(), variantReport.consequence());
        assertEquals(COSMIC_ID, variantReport.cosmicID());
        assertEquals(TOTAL_READ_COUNT, variantReport.totalReadCount());
        assertEquals(ALLELE_READ_COUNT, variantReport.alleleReadCount());
    }
}
