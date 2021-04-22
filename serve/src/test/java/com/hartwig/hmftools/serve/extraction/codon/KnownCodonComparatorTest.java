package com.hartwig.hmftools.serve.extraction.codon;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.junit.Test;

public class KnownCodonComparatorTest {

    @Test
    public void canSortKnownCodons() {
        KnownCodon codon1 = ImmutableKnownCodon.builder()
                .annotation(ImmutableCodonAnnotation.builder()
                        .chromosome("1")
                        .start(10)
                        .end(11)
                        .gene("gene x")
                        .mutationType(MutationTypeFilter.ANY)
                        .codonIndex(1)
                        .transcript("transcript x")
                        .build())
                .build();

        KnownCodon codon2 = ImmutableKnownCodon.builder()
                .annotation(ImmutableCodonAnnotation.builder().from(codon1.annotation()).codonIndex(2).build())
                .build();

        Set<KnownCodon> sortedExons = Sets.newTreeSet(new KnownCodonComparator());
        sortedExons.add(codon2);
        sortedExons.add(codon1);

        assertEquals(codon1, sortedExons.iterator().next());
    }

}