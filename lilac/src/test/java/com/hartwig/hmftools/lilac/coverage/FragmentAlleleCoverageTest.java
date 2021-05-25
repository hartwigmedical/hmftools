package com.hartwig.hmftools.lilac.coverage;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class FragmentAlleleCoverageTest
{
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
                createFragment("01"), Lists.newArrayList(allele1), Lists.newArrayList(), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("02"), Lists.newArrayList(allele2), Lists.newArrayList(), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("03"), Lists.newArrayList(allele3), Lists.newArrayList(), Lists.newArrayList()));

        fragmentAlleles.add(new FragmentAlleles(
                createFragment("04"), Lists.newArrayList(allele2), Lists.newArrayList(allele1), Lists.newArrayList(allele3)));

        // wild-only is ignored
        fragmentAlleles.add(new FragmentAlleles(
                createFragment("04"), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(allele1, allele2, allele3)));

        FragmentAlleleMatrix matrix = new FragmentAlleleMatrix(fragmentAlleles, alleles);

        List<HlaAlleleCoverage> coverages = matrix.create(complex);
        assertEquals(3, coverages.size());
        assertEquals(1.33, coverages.get(0).TotalCoverage, 0.01);
        assertEquals(1.33, coverages.get(1).TotalCoverage, 0.01);
        assertEquals(1.33, coverages.get(2).TotalCoverage, 0.01);
        assertEquals(1, coverages.get(0).UniqueCoverage, 0.01);
        assertEquals(1, coverages.get(1).UniqueCoverage, 0.01);
        assertEquals(1, coverages.get(2).UniqueCoverage, 0.01);
        assertEquals(0, coverages.get(0).WildCoverage, 0.01);
        assertEquals(0, coverages.get(1).WildCoverage, 0.01);
        assertEquals(0.33, coverages.get(2).WildCoverage, 0.01);
    }

    private AminoAcidFragment createFragment(final String id)
    {
        return new AminoAcidFragment(
                id, "", Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(),
                Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());
    }
}
