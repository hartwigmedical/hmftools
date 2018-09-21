package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import com.google.common.collect.Lists;

import org.junit.Test;

public class CanonicalAnnotationTest {

    private static final String FIELD_SEPARATOR = "\\|";

    @Test
    public void favourDriverCatalogGenesRegardlessOfCodingEffect() {

        final String first  = "C|frameshift_variant|HIGH|ATP1A1|ENSG00000163399|transcript|ENST00000295598|protein_coding|2/23|c.74delA|p.Lys25fs|326/3654|74/3072|25/1023||INFO_REALIGN_3_PRIME";
        final String second = "C|frameshift_variant|HIGH|AL136376.1|ENSG00000269279|transcript|ENST00000598661|protein_coding|1/2|c.26delT|p.Leu9fs|26/342|26/342|9/113||WARNING_TRANSCRIPT_NO_STOP_CODON";
        final String third  = "C|synonymous_variant|HIGH|ATP1A1|ENSG00000163399|transcript|ENST00000537345|protein_coding|2/23|c.74delA|p.Lys25fs|437/3763|74/3072|25/1023||INFO_REALIGN_3_PRIME";

        final SnpEffAnnotation unknownGene = SnpEffAnnotationFactory.fromParts(first.split(FIELD_SEPARATOR));
        final SnpEffAnnotation knownButNotDriverGene = SnpEffAnnotationFactory.fromParts(second.split(FIELD_SEPARATOR));
        final SnpEffAnnotation synonymousDriverGene = SnpEffAnnotationFactory.fromParts(third.split(FIELD_SEPARATOR));

        final CanonicalAnnotation victim = new CanonicalAnnotation();
        assertEquals(Optional.empty(), victim.canonicalSnpEffAnnotation(Lists.newArrayList(unknownGene)));
        assertEquals(knownButNotDriverGene, victim.canonicalSnpEffAnnotation(Lists.newArrayList(unknownGene, knownButNotDriverGene)).get());
        assertEquals(synonymousDriverGene, victim.canonicalSnpEffAnnotation(Lists.newArrayList(unknownGene, knownButNotDriverGene, synonymousDriverGene)).get());
    }

}
