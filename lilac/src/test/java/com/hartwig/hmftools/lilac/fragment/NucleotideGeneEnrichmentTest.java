package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacUtils.namesMatch;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragment.expandIndices;

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
        final List<Integer> A_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 364, 365);
        final List<Integer> B_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 337, 348, 362);
        final List<Integer> C_EXON_BOUNDARIES = Lists.newArrayList(24, 114, 206, 298, 338, 349, 365, 366);

        NucleotideGeneEnrichment enricher = new NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);

        List<Integer> indices = Lists.newArrayList();
        indices.add(337);
        indices.add(348);
        assertGene(enricher, Sets.newHashSet("HLA-A", "HLA-B"), "HLA-A", indices);

        assertGene(enricher, Sets.newHashSet("HLA-A", "HLA-B"), "HLA-B", indices);
        assertGene(enricher, Sets.newHashSet("HLA-C"), "HLA-C", indices);

        indices.add(362);
        assertGene(enricher, Sets.newHashSet("HLA-A"), "HLA-A", indices);
        assertGene(enricher, Sets.newHashSet("HLA-B"), "HLA-B", indices);
    }

    private void assertGene(
            NucleotideGeneEnrichment enricher,
            Set<String> expectedGenes, String alignedGene, List<Integer> aminoAcideIndices)
    {
        NucleotideFragment victim = create(alignedGene, expandIndices(aminoAcideIndices));
        NucleotideFragment result = enricher.enrich(victim);
        assertTrue(namesMatch(result.getGenes(), expectedGenes));
    }

    private NucleotideFragment create(String gene, List<Integer> indices)
    {
        List<Integer> qualities = Lists.newArrayList();
        qualities.add(0);
        return new NucleotideFragment("id", "", Sets.newHashSet(gene), indices, qualities, Lists.newArrayList("G"));
    }

}
