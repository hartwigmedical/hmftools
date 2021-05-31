package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.WILD_ONLY;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildLoci;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildTargetSequences;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.junit.Test;

public class FragmentAllelesTest
{
    @Test
    public void testFragmentAlleleAssignment()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaSequenceLoci seq1 = new HlaSequenceLoci(allele1, Lists.newArrayList("A", "B", "C", "D"));
        List<HlaSequenceLoci> sequences = Lists.newArrayList(seq1);

        String readInfo = "";
        List<String> emptyNucs = Lists.newArrayList();
        List<Integer> emptyQuals = Lists.newArrayList();
        List<Integer> emptyLoci = Lists.newArrayList();
        Set<String> aGenes = Sets.newHashSet(HLA_A);
        Set<String> bGenes = Sets.newHashSet(HLA_B);

        final Map<String,Map<Integer,List<String>>> geneHetLociMap = Maps.newHashMap();
        Map<Integer,List<String>> geneALociMap = Maps.newLinkedHashMap();
        geneALociMap.put(0, Lists.newArrayList("A"));
        geneALociMap.put(1, Lists.newArrayList("B"));
        geneALociMap.put(2, Lists.newArrayList("C"));
        geneALociMap.put(3, Lists.newArrayList("D"));
        geneHetLociMap.put(GENE_A, geneALociMap);

        Map<String,List<Integer>> refNucleotideHetLoci = Maps.newHashMap();
        List<HlaSequenceLoci> candidateNucSequences = Lists.newArrayList();
        List<Set<String>> refNucleotides = Lists.newArrayList();

        // basic match
        Fragment frag1 = new Fragment("01", readInfo, aGenes, emptyLoci, emptyQuals, emptyNucs);
        frag1.setAminoAcids(Lists.newArrayList(0, 2, 3), Lists.newArrayList("A", "C", "D"));
        List<Fragment> fragments = Lists.newArrayList(frag1);

        FragmentAlleleMapper mapper = new FragmentAlleleMapper(geneHetLociMap, refNucleotideHetLoci, refNucleotides);
        List<FragmentAlleles> fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertEquals(1, fragAlleles.size());
        assertTrue(fragAlleles.get(0).getFull().contains(allele1));

        // check non-het locations aren't checked
        Map<Integer,List<String>> geneBLociMap = Maps.newLinkedHashMap();
        geneBLociMap.put(1, Lists.newArrayList("B"));
        geneBLociMap.put(3, Lists.newArrayList("G"));
        geneHetLociMap.put(GENE_B, geneBLociMap);

        HlaAllele allele2 = HlaAllele.fromString("B*01:01");
        HlaSequenceLoci seq2 = new HlaSequenceLoci(allele2, Lists.newArrayList("A", "B", "C", "D"));
        sequences = Lists.newArrayList(seq2);

        Fragment frag2 = new Fragment("02", readInfo, bGenes, emptyLoci, emptyQuals, emptyNucs);
        frag2.setAminoAcids(Lists.newArrayList(0, 1, 2, 3), Lists.newArrayList("B", "B", "D", "D"));
        fragments = Lists.newArrayList(frag2);

        fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertEquals(1, fragAlleles.size());
        assertTrue(fragAlleles.get(0).getFull().contains(allele2));

        // wildcard match
        HlaAllele allele3 = HlaAllele.fromString("A*01:22");
        HlaSequenceLoci seq3 = new HlaSequenceLoci(allele3, Lists.newArrayList("*", "B", "C", "*"));
        sequences = Lists.newArrayList(seq3);

        Fragment frag3 = new Fragment("03", readInfo, aGenes, emptyLoci, emptyQuals, emptyNucs);
        frag3.setAminoAcids(Lists.newArrayList(0, 1, 2, 3), Lists.newArrayList("A", "B", "C", "D"));
        fragments = Lists.newArrayList(frag3);

        fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertTrue(fragAlleles.isEmpty());
        assertEquals(WILD_ONLY, frag3.scope());

        // when a wild is coupled with a full match (non-wild), it is retained
        sequences = Lists.newArrayList(seq1, seq3);
        fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertEquals(1, fragAlleles.size());
        assertTrue(fragAlleles.get(0).getFull().contains(allele1));
        assertTrue(fragAlleles.get(0).getWild().contains(allele3));
    }

    @Test
    public void testFragmentAlleleLowBaseQuals()
    {
        // recover low-qual bases to establish a match
        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaSequenceLoci seq1 = new HlaSequenceLoci(allele1, Lists.newArrayList("A", "P", "C", "D"));
        List<HlaSequenceLoci> sequences = Lists.newArrayList(seq1);

        String readInfo = "";
        Set<String> aGenes = Sets.newHashSet(HLA_A);

        final Map<String,Map<Integer,List<String>>> geneHetLociMap = Maps.newHashMap();
        Map<Integer,List<String>> geneALociMap = Maps.newLinkedHashMap();
        geneALociMap.put(0, Lists.newArrayList("A"));
        geneALociMap.put(1, Lists.newArrayList("P"));
        geneALociMap.put(2, Lists.newArrayList("C"));
        geneALociMap.put(3, Lists.newArrayList("D"));
        geneHetLociMap.put(GENE_A, geneALociMap);

        Map<String,List<Integer>> refNucleotideHetLoci = Maps.newHashMap();
        List<HlaSequenceLoci> candidateNucSequences = Lists.newArrayList();
        List<Set<String>> refNucleotides = Lists.newArrayList();

        // GCT -> A, CCC -> P, TGC -> C, GAC -> D
        String sequence = "GCTCCCTGCGAC";
        List<Integer> nucLoci = buildLoci(sequence);
        List<String> nucleotides = buildTargetSequences(sequence, nucLoci);
        List<Integer> nucQuals = nucLoci.stream().map(x -> 30).collect(Collectors.toList());

        // mark 1st, 2nd and 4th AAs with low-qual bases
        nucQuals.set(1, 2);
        nucQuals.set(4, 2);
        nucQuals.set(11, 2);
        Fragment frag1 = new Fragment("01", readInfo, aGenes, nucLoci, nucQuals, nucleotides);

        assertTrue(frag1.validate());

        frag1.qualityFilter(30);
        frag1.buildAminoAcids();
        assertEquals(9, frag1.getNucleotideLoci().size());
        assertEquals(1, frag1.getAminoAcidLoci().size());

        List<Fragment> fragments = Lists.newArrayList(frag1);

        FragmentAlleleMapper mapper = new FragmentAlleleMapper(geneHetLociMap, refNucleotideHetLoci, refNucleotides);
        List<FragmentAlleles> fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertEquals(1, fragAlleles.size());
        assertTrue(fragAlleles.get(0).getFull().contains(allele1));

        // again but with the low-qual bases not forming the correct amino acid, but still accepted as a match
        frag1.getRawNucleotides().set(1, "A");

        fragAlleles = mapper.createFragmentAlleles(fragments, sequences, candidateNucSequences);

        assertEquals(1, fragAlleles.size());
        assertTrue(fragAlleles.get(0).getFull().contains(allele1));
    }

    @Test
    public void testFragmentAlleleCoverage()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaAllele allele2 = HlaAllele.fromString("B*01:01");
        HlaAllele allele3 = HlaAllele.fromString("C*01:01");

        List<HlaAllele> alleles = Lists.newArrayList(allele1, allele2, allele3);
        HlaComplex complex = new HlaComplex(alleles);

        List<FragmentAlleles> fragmentAlleles = Lists.newArrayList();

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("01"), Lists.newArrayList(allele1), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("02"), Lists.newArrayList(allele2), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("03"), Lists.newArrayList(allele3), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("04"), Lists.newArrayList(allele2), Lists.newArrayList(allele1, allele3)));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("05"), Lists.newArrayList(), Lists.newArrayList(allele1, allele2, allele3)));

        FragmentAlleleMatrix matrix = new FragmentAlleleMatrix(fragmentAlleles, alleles);

        List<HlaAlleleCoverage> coverages = matrix.create(complex);
        assertEquals(3, coverages.size());
        assertEquals(1.67, coverages.get(0).TotalCoverage, 0.01);
        assertEquals(1.67, coverages.get(1).TotalCoverage, 0.01);
        assertEquals(1.67, coverages.get(2).TotalCoverage, 0.01);
        assertEquals(1, coverages.get(0).UniqueCoverage, 0.01);
        assertEquals(1, coverages.get(1).UniqueCoverage, 0.01);
        assertEquals(1, coverages.get(2).UniqueCoverage, 0.01);
        assertEquals(0.67, coverages.get(0).WildCoverage, 0.01);
        assertEquals(0.33, coverages.get(1).WildCoverage, 0.01);
        assertEquals(0.67, coverages.get(2).WildCoverage, 0.01);
    }

    private Fragment createFragment(final String id)
    {
        return new Fragment(
                id, "", Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());
    }

    public static class HlaComplexTest
    {
        HlaAllele a11 = createA(1, 1);
        HlaAllele a12 = createA(1, 2);
        HlaAllele a13 = createA(1, 3);
        HlaAllele a21 = createA(2, 1);
        HlaAllele a22 = createA(2, 2);
        HlaAllele a23 = createA(2, 3);
        HlaAllele a31 = createA(3, 1);

        List<HlaAllele> all = Lists.newArrayList(a11, a12, a13, a21, a22, a23, a31);

        HlaAllele b11 = createB(1, 1);
        HlaAllele b12 = createB(1, 2);
        HlaAllele b13 = createB(1, 3);

        @Test
        public void testTwoConfirmedProtein()
        {
            List<HlaAllele> confirmedProtein = Lists.newArrayList(a11, a12);
            List<HlaAllele> confirmedGroup = confirmedProtein.stream().map(x -> x.asAlleleGroup()).collect(Collectors.toList());
            // distinct

            List<HlaComplex> result = HlaComplexBuilder.buildComplexesByGene(GENE_A, confirmedGroup, confirmedProtein, all);
            assertEquals(confirmedProtein, result.get(0).getAlleles());
        }

        /* print only
        @Test
        public void testOneConfirmedGroup()
        {
            List<HlaAllele> confirmedGroup = Lists.newArrayList(a11.asAlleleGroup());

            List<HlaComplex> result = HlaComplex.gene("A", confirmedGroup, Collections.emptyList(), all);

            for(hlaComplex in result)
            {
                println(hlaComplex)
            }

        }

        @Test
        public void testNoConfirmed()
        {
            val result = HlaComplex.gene("A", Collections.emptyList(), Collections.emptyList(), all)
            for(hlaComplex in result)
            {
                println(hlaComplex)
            }

        }

        @Test
        public void testCombineComplexes()
        {
            val complexA1 = HlaComplex(listOf(a11, a12))
            val complexA2 = HlaComplex(listOf(a11, a13))
            val complexA3 = HlaComplex(listOf(a12, a13))
            val complexA = listOf(complexA1, complexA2, complexA3)

            val complexB1 = HlaComplex(listOf(b11, b12))
            val complexB2 = HlaComplex(listOf(b11, b13))
            val complexB = listOf(complexB1, complexB2)

            val result = HlaComplex.combineComplexes(complexA, complexB)
            for(hlaComplex in result)
            {
                println(hlaComplex.alleles)
            }

        }

         */

        private HlaAllele createA(int group, int protein)
        {
            return new HlaAllele(GENE_A, String.valueOf(group), String.valueOf(protein), "", "", null, null);
        }

        private HlaAllele createB(int group, int protein)
        {
            return new HlaAllele(GENE_B, String.valueOf(group), String.valueOf(protein), "", "", null, null);
        }

    }
}
