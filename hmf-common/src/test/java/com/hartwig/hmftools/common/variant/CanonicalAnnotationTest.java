package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.cosmic.CosmicAnnotation;
import com.hartwig.hmftools.common.variant.cosmic.ImmutableCosmicAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.ImmutableSnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CanonicalAnnotationTest {

    private static final String CDKN2A = "ENST00000498124";
    private static final String CDKN2Ap14ARF = "ENST00000361570";

    @Test
    public void testFavourCDKN2ACosmicAnnotation() {
        final CosmicAnnotation p16 = createCosmicAnnotation("CDKN2A", CDKN2A);
        final CosmicAnnotation p14 = createCosmicAnnotation("CDKN2Ap14ARF", CDKN2Ap14ARF);

        final List<CosmicAnnotation> all = Lists.newArrayList(p14, p16);
        Collections.shuffle(all);

        final CanonicalAnnotation victim = new CanonicalAnnotation();
        assertEquals(p16, victim.canonicalCosmicAnnotation(all).get());
        assertEquals(p14, victim.canonicalCosmicAnnotation(Lists.newArrayList(p14)).get());
    }

    @Test
    public void testFavourCDKN2ASnpEffAnnotation() {
        final SnpEffAnnotation p16 =
                createSnpEffAnnotation("CDKN2A", CDKN2A, Lists.newArrayList(VariantConsequence.MISSENSE_VARIANT));
        final SnpEffAnnotation p14 =
                createSnpEffAnnotation("CDKN2Ap14ARF", CDKN2Ap14ARF, Lists.newArrayList(VariantConsequence.MISSENSE_VARIANT));

        final List<SnpEffAnnotation> all = Lists.newArrayList(p14, p16);
        Collections.shuffle(all);

        final CanonicalAnnotation victim = new CanonicalAnnotation();
        assertEquals(p16, victim.canonicalSnpEffAnnotation(all).get());
        assertEquals(p14, victim.canonicalSnpEffAnnotation(Lists.newArrayList(p14)).get());
    }

    @Test
    public void favourDriverCatalogGenes() {
        final TranscriptAnnotation nonCanonicalDriverGene = createCosmicAnnotation("ATP1A1", "ENST00000295598");
        final TranscriptAnnotation noDriverGene = createCosmicAnnotation("AL136376.1", "ENST00000598661");
        final TranscriptAnnotation canonicalDriverGene = createCosmicAnnotation("ATP1A1", "ENST00000537345");

        final CanonicalAnnotation victim = new CanonicalAnnotation();
        assertEquals(Optional.empty(), victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene)));

        Optional<TranscriptAnnotation> annotationSecond =
                victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene, noDriverGene));
        assertTrue(annotationSecond.isPresent());
        assertEquals(noDriverGene, annotationSecond.get());

        Optional<TranscriptAnnotation> annotationThird =
                victim.pickCanonicalFavourDriverGene(Lists.newArrayList(nonCanonicalDriverGene, noDriverGene, canonicalDriverGene));
        assertTrue(annotationThird.isPresent());
        assertEquals(canonicalDriverGene, annotationThird.get());
    }

    @NotNull
    private static CosmicAnnotation createCosmicAnnotation(@NotNull final String gene, @NotNull final String transcript) {
        return ImmutableCosmicAnnotation.builder()
                .gene(gene)
                .transcript(transcript)
                .count(0)
                .id("1")
                .hgvsCoding("")
                .hgvsProtein("")
                .build();
    }

    @NotNull
    private static SnpEffAnnotation createSnpEffAnnotation(@NotNull final String gene, @NotNull final String transcript,
            @NotNull final List<VariantConsequence> consequences) {
        return ImmutableSnpEffAnnotation.builder()
                .allele("")
                .effects("")
                .consequences(consequences)
                .severity("")
                .gene(gene)
                .geneID(gene)
                .featureType("transcript")
                .featureID(transcript)
                .transcriptBioType("")
                .rank("")
                .hgvsCoding("")
                .hgvsProtein("")
                .cDNAPosAndLength("")
                .cdsPosAndLength("")
                .aaPosAndLength("")
                .distance("")
                .addition("")
                .build();
    }
}
