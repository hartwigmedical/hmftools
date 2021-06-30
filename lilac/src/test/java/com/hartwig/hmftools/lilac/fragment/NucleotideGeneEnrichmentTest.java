package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacUtils.namesMatch;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;

import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class NucleotideGeneEnrichmentTest
{
    @Test
    public void testGeneEnrichment()
    {
        NucleotideGeneEnrichment enricher = new NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);

        List<Integer> indices = Lists.newArrayList();
        indices.add(337);
        indices.add(348);
        assertGene(enricher, Sets.newHashSet(HLA_A, HLA_B), HLA_A, indices);

        assertGene(enricher, Sets.newHashSet(HLA_A, HLA_B), HLA_B, indices);
        assertGene(enricher, Sets.newHashSet(HLA_C), HLA_C, indices);

        indices.add(365);
        assertGene(enricher, Sets.newHashSet("HLA-A"), "HLA-A", indices);
        assertGene(enricher, Sets.newHashSet("HLA-B"), "HLA-B", indices);
    }

    private void assertGene(
            NucleotideGeneEnrichment enricher,
            Set<String> expectedGenes, String alignedGene, List<Integer> aminoAcideIndices)
    {
        Fragment fragment = create(alignedGene, expandIndices(aminoAcideIndices));
        Fragment result = enricher.enrich(fragment);
        assertTrue(namesMatch(result.getGenes(), expectedGenes));
    }

    private Fragment create(final String gene, final List<Integer> indices)
    {
        List<Integer> qualities = Lists.newArrayList();
        qualities.add(0);
        return new Fragment("id", "", gene, Sets.newHashSet(gene), indices, qualities, Lists.newArrayList("G"));
    }

}
