package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class EnrichedSomaticVariantFactoryTest {

    @Test
    public void revertToNormalCosmicIDIfCanonicalIsNoMatch() {
        final String line =
                "14\t105246551\trs121434592;COSM33765\tC\tT\t.\tPASS\tCOSM2ENST=COSM33765|AKT1_ENST00000349310|c.49G>A|p.E17K|519\tGT:AD:DP\t0/1:120,71:204";

        final String sample = "sample";
        final VariantContext context = VariantContextFromString.decode(sample, line);
        final Optional<SomaticVariant> optVariant = SomaticVariantFactory.passOnlyInstance().createVariant(sample, context);
        assert optVariant.isPresent();
        final SomaticVariant variant = optVariant.get();

        final ImmutableEnrichedSomaticVariant.Builder builder = SomaticVariantTestBuilderFactory.createEnriched();
        final Map<String, String> geneToTranscriptMap = Maps.newHashMap();
        geneToTranscriptMap.put("AKT1", "ENST0123456");
        final TranscriptAnnotationSelector selector = new TranscriptAnnotationSelector(geneToTranscriptMap);

        EnrichedSomaticVariantFactory.addCanonicalCosmicID(builder, variant, selector);

        final EnrichedSomaticVariant enrichedSomaticVariant = builder.build();

        assertEquals("COSM33765", enrichedSomaticVariant.canonicalCosmicID());
    }

    @Test
    public void testCanonicalCodingEffectUsesTranscriptAnnotation() {
        final String line =
                "11\t133715264\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|sequence_feature|MODERATE|SPATA19|ENSG00000166118|modified-residue:Phosphoserine|ENST00000299140|protein_coding|1/7|c.78G>A||||||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000299140|protein_coding|1/7|c.78G>A|p.Ser26Ser|133/861|78/504|26/167||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000532889|protein_coding|1/7|c.78G>A|p.Ser26Ser|170/653|78/504|26/167||\tGT:AD:DP\t0/1:57,49:108";
        final String sample = "sample";
        final VariantContext context = VariantContextFromString.decode(sample, line);
        final Optional<SomaticVariant> optVariant = SomaticVariantFactory.passOnlyInstance().createVariant(sample, context);
        assert optVariant.isPresent();
        final SomaticVariant variant = optVariant.get();

        final ImmutableEnrichedSomaticVariant.Builder builder = SomaticVariantTestBuilderFactory.createEnriched();
        final Map<String, String> geneToTranscriptMap = Maps.newHashMap();
        geneToTranscriptMap.put("SPATA19", "ENST00000299140");

        final TranscriptAnnotationSelector selector = new TranscriptAnnotationSelector(geneToTranscriptMap);
        EnrichedSomaticVariantFactory.addCanonicalEffect(builder, variant, selector);

        final EnrichedSomaticVariant enrichedSomaticVariant = builder.build();

        assertEquals(CodingEffect.SYNONYMOUS, enrichedSomaticVariant.canonicalCodingEffect());
    }

}
