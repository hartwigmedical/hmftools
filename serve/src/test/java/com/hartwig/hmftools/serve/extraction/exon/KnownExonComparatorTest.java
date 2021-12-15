package com.hartwig.hmftools.serve.extraction.exon;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.junit.Test;

public class KnownExonComparatorTest {

    @Test
    public void canSortKnownExons() {
        KnownExon exon1 = ImmutableKnownExon.builder()
                .annotation(ImmutableExonAnnotation.builder()
                        .chromosome("1")
                        .start(10)
                        .end(11)
                        .gene("gene x")
                        .mutationType(MutationTypeFilter.ANY)
                        .rank(1)
                        .transcript("transcript x")
                        .build())
                .build();

        KnownExon exon2 = ImmutableKnownExon.builder()
                .annotation(ImmutableExonAnnotation.builder().from(exon1.annotation()).rank(2).build())
                .build();

        Set<KnownExon> sortedExons = Sets.newTreeSet(new KnownExonComparator());
        sortedExons.add(exon2);
        sortedExons.add(exon1);

        assertEquals(exon1, sortedExons.iterator().next());
    }
}