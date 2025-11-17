package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.LilacUtils.namesMatch;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.app.LilacAppTest.buildGeneCache;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.expandIndices;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createReadRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.ArgumentMatchers.anyInt;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.NavigableMap;
import java.util.Set;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConstants;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.util.ThrowOnUnstubbed;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

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

        NucleotideGeneEnrichment enricher = NucleotideGeneEnrichment.create(GENE_EXON_BOUNDARIES);

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

    @Test
    public void testApplyQualityFilterSingleNucleotideFilteredOut()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        List<Nucleotide> fragmentNucleotides = Lists.newArrayList(
                Nucleotide.create(0, (byte) 40, "A"),
                Nucleotide.create(1, (byte) 40, "T"));
        final NavigableMap<Integer, Nucleotide> fragmentNucleotidesByLoci = Maps.newTreeMap();
        fragmentNucleotides.forEach(x -> fragmentNucleotidesByLoci.put(x.locus(), x));

        SAMRecord samRecord = mock(SAMRecord.class);
        Read read = mock(Read.class);
        doReturn(samRecord).when(read).bamRecord();
        List<Read> reads = Lists.newArrayList(read);

        Fragment fragment = mock(Fragment.class, new ThrowOnUnstubbed());
        doReturn(fragmentNucleotidesByLoci).when(fragment).nucleotidesByLoci();
        doReturn(reads).when(fragment).reads();
        doReturn(HLA_A).when(fragment).readGene();
        doReturn(Sets.newHashSet(HLA_A)).when(fragment).genes();

        List<String> minEvidenceSequences = Lists.newArrayList("A");
        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(minEvidenceSequences).when(pooledCounts).getMinEvidenceSequences(anyInt());

        Multiset<String> localNucleotides = HashMultiset.create();
        localNucleotides.setCount("A", DEFAULT_MIN_DEPTH_FILTER);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localNucleotides).when(localCounts).get(anyInt());

        Fragment result = NucleotideFragmentQualEnrichment.applyQualityFilter(fragment, pooledCounts, pooledCounts, localCounts);
        NavigableMap<Integer, Nucleotide> actualNucleotides = result.nucleotidesByLoci();
        NavigableMap<Integer, Nucleotide> expectedNucleotides = Maps.newTreeMap();
        expectedNucleotides.put(0, fragmentNucleotides.get(0));

        assertEquals(expectedNucleotides, actualNucleotides);
    }

    @Test
    public void testApplyQualityFilterSingleNucleotideSavedDueToLowDepth()
    {
        LilacConstants.MIN_DEPTH_FILTER = DEFAULT_MIN_DEPTH_FILTER;

        List<Nucleotide> fragmentNucleotides = Lists.newArrayList(
                Nucleotide.create(0, (byte) 40, "A"),
                Nucleotide.create(1, (byte) 40, "T"));
        final NavigableMap<Integer, Nucleotide> fragmentNucleotidesByLoci = Maps.newTreeMap();
        fragmentNucleotides.forEach(x -> fragmentNucleotidesByLoci.put(x.locus(), x));

        SAMRecord samRecord = mock(SAMRecord.class);
        Read read = mock(Read.class);
        doReturn(samRecord).when(read).bamRecord();
        List<Read> reads = Lists.newArrayList(read);

        Fragment fragment = mock(Fragment.class, new ThrowOnUnstubbed());
        doReturn(fragmentNucleotidesByLoci).when(fragment).nucleotidesByLoci();
        doReturn(reads).when(fragment).reads();
        doReturn(HLA_A).when(fragment).readGene();
        doReturn(Sets.newHashSet(HLA_A)).when(fragment).genes();

        List<String> minEvidenceSequences = Lists.newArrayList("A");
        SequenceCount pooledCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(minEvidenceSequences).when(pooledCounts).getMinEvidenceSequences(anyInt());

        Multiset<String> localHighNucleotides = HashMultiset.create();
        localHighNucleotides.setCount("A", DEFAULT_MIN_DEPTH_FILTER);
        Multiset<String> localLowNucleotides = HashMultiset.create();
        localHighNucleotides.setCount("A", DEFAULT_MIN_DEPTH_FILTER - 1);
        SequenceCount localCounts = mock(SequenceCount.class, new ThrowOnUnstubbed());
        doReturn(localHighNucleotides).when(localCounts).get(0);
        doReturn(localLowNucleotides).when(localCounts).get(1);

        Fragment result = NucleotideFragmentQualEnrichment.applyQualityFilter(fragment, pooledCounts, pooledCounts, localCounts);
        NavigableMap<Integer, Nucleotide> actualNucleotides = result.nucleotidesByLoci();

        assertEquals(fragmentNucleotidesByLoci, actualNucleotides);
    }

    private static void assertGene(
            final NucleotideGeneEnrichment enricher,
            final Set<HlaGene> expectedGenes, final HlaGene alignedGene, final Iterable<Integer> aminoAcideIndices)
    {
        Fragment fragment = create(alignedGene, expandIndices(aminoAcideIndices));

        enricher.checkAdditionalGenes(fragment);
        assertTrue(namesMatch(fragment.genes(), expectedGenes));
    }

    private static Fragment create(final HlaGene gene, final List<Integer> indices)
    {
        List<Byte> qualities = Lists.newArrayListWithCapacity(indices.size());
        List<String> nucleotides = Lists.newArrayListWithCapacity(indices.size());

        for(int i = 0; i < indices.size(); ++i)
        {
            qualities.add((byte) 0);
            nucleotides.add("G");
        }

        return Fragment.createFromQuals(createReadRecord("01"), gene, Sets.newHashSet(gene), indices, qualities, nucleotides);
    }
}
