package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.junit.Test;

public class FragmentAllelesTest
{
    @Test
    public void testFragmentAlleleCreation()
    {
        /*
        final List<Fragment> refCoverageFragments, final List<Integer> refAminoAcidHetLoci,
        final List<HlaSequenceLoci> candidateAminoAcidSequences, final List<Set<String>> refAminoAcids,
        final Map<String,List<Integer>> refNucleotideHetLoci, final List<HlaSequenceLoci> candidateNucleotideSequences,
        final List<Set<String>> refNucleotides)



        List<FragmentAlleles> fragAlleles = FragmentAlleles.createFragmentAlleles(
        final List<Fragment> refCoverageFragments, final List<Integer> refAminoAcidHetLoci,
        final List<HlaSequenceLoci> candidateAminoAcidSequences, final List<Set<String>> refAminoAcids,
        final Map<String,List<Integer>> refNucleotideHetLoci, final List<HlaSequenceLoci> candidateNucleotideSequences,
        final List<Set<String>> refNucleotides)
*/
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
