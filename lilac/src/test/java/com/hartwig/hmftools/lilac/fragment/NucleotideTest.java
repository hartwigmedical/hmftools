package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacUtils.namesMatch;
import static com.hartwig.hmftools.lilac.ReferenceData.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.ReferenceData.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.ReferenceData.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.app.LilacAppTest.buildGeneCache;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createReadRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class NucleotideTest
{
    @Test
    public void testCreateNucleotidesFromAminoAcid()
    {
        assertEquals(Lists.newArrayList("T", "A", "A"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("X"));
        assertEquals(Lists.newArrayList(".", ".", "."), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("."));
        assertEquals(Lists.newArrayList("A", "G", "CTAA"), NucleotideFragmentFactory.createNucleotidesFromAminoAcid("SX"));
    }

    @Test
    public void testGeneEnrichment()
    {
        buildGeneCache();

        NucleotideGeneEnrichment enricher = new NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);

        List<Integer> indices = Lists.newArrayList();
        indices.add(337);
        indices.add(348);
        assertGene(enricher, Sets.newHashSet(HLA_A, HLA_B), HLA_A, indices);

        assertGene(enricher, Sets.newHashSet(HLA_A, HLA_B), HLA_B, indices);
        assertGene(enricher, Sets.newHashSet(HLA_C), HLA_C, indices);

        indices.add(365);
        assertGene(enricher, Sets.newHashSet(HLA_A), HLA_A, indices);
        assertGene(enricher, Sets.newHashSet(HLA_B), HLA_B, indices);
    }

    private static void assertGene(
            final NucleotideGeneEnrichment enricher,
            final Set<String> expectedGenes, final String alignedGene, final List<Integer> aminoAcideIndices)
    {
        Fragment fragment = create(alignedGene, expandIndices(aminoAcideIndices));

        enricher.checkAdditionalGenes(fragment);
        assertTrue(namesMatch(fragment.genes(), expectedGenes));
    }

    private static Fragment create(final String gene, final List<Integer> indices)
    {
        List<Byte> qualities = Lists.newArrayListWithCapacity(indices.size());
        List<String> nucleotides = Lists.newArrayListWithCapacity(indices.size());

        for(int i = 0; i < indices.size(); ++i)
        {
            qualities.add((byte) 0);
            nucleotides.add("G");
        }

        return new Fragment(createReadRecord("01"), gene, Sets.newHashSet(gene), indices, qualities, nucleotides);
    }

}
