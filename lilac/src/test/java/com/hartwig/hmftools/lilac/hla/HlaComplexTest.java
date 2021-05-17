package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.coverage.HlaComplex;

import org.junit.Test;

public class HlaComplexTest
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

        List<HlaComplex> result = HlaComplex.gene(GENE_A, confirmedGroup, confirmedProtein, all);
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
