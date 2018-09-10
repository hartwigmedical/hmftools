package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.snpeff.AnnotationTestFactory.createVariantAnnotationBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.ImmutableSnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

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

    private static final String CODING_IMPACT = "c.RtoA";
    private static final String PROTEIN_IMPACT = "p.RtoA";
    private static final String VARIANT_DETAILS = CODING_IMPACT + " (" + PROTEIN_IMPACT + ")";

    @Test
    public void worksAsExpected() {
        final ConsequenceDeterminer determiner = new ConsequenceDeterminer(Sets.newHashSet(TRANSCRIPT));

        final VariantConsequence rightConsequence = VariantConsequence.MISSENSE_VARIANT;
        final VariantConsequence wrongConsequence = VariantConsequence.OTHER;

        final ImmutableSnpEffAnnotation.Builder annotationBuilder = createVariantAnnotationBuilder().featureID(TRANSCRIPT).
                featureType(ConsequenceDeterminer.FEATURE_TYPE_TRANSCRIPT).gene(GENE).hgvsCoding(CODING_IMPACT).
                hgvsProtein(PROTEIN_IMPACT);
        final SnpEffAnnotation rightAnnotation = annotationBuilder.consequences(Lists.newArrayList(rightConsequence)).build();
        final SnpEffAnnotation wrongAnnotation = annotationBuilder.consequences(Lists.newArrayList(wrongConsequence)).build();

        final ImmutableEnrichedSomaticVariant.Builder variantBuilder = SomaticVariantTestBuilderFactory.createEnriched().
                chromosome(CHROMOSOME).ref(REF).alt(ALT).canonicalCosmicID(COSMIC_ID).position(POSITION).
                totalReadCount(TOTAL_READ_COUNT).alleleReadCount(ALLELE_READ_COUNT);

        final EnrichedSomaticVariant rightVariant = variantBuilder.snpEffAnnotations(Lists.newArrayList(rightAnnotation)).build();
        final EnrichedSomaticVariant wrongVariant = variantBuilder.snpEffAnnotations(Lists.newArrayList(wrongAnnotation)).build();

        final List<VariantReport> variantReports = determiner.run(Lists.newArrayList(rightVariant, wrongVariant));
        assertEquals(1, variantReports.size());

        final VariantReport variantReport = variantReports.get(0);
        assertEquals(GENE, variantReport.gene());
        assertEquals(VARIANT_DETAILS, variantReport.variantDetails());
        assertEquals(COSMIC_ID, variantReport.cosmicID());
        assertEquals(TOTAL_READ_COUNT, variantReport.totalReadCount());
        assertEquals(ALLELE_READ_COUNT, variantReport.alleleReadCount());
    }
}
