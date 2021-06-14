package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class CoverageTest
{
    @Test
    public void testAlleleCoverageExpansion()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaAllele allele2 = HlaAllele.fromString("A*02:01");
        HlaAllele allele3 = HlaAllele.fromString("B*01:01");
        HlaAllele allele4 = HlaAllele.fromString("B*02:01");
        HlaAllele allele5 = HlaAllele.fromString("C*01:01");
        HlaAllele allele6 = HlaAllele.fromString("C*02:01");

        // standard 6-allele solution
        HlaAlleleCoverage cov1 = new HlaAlleleCoverage(allele1, 10,20, 0);
        HlaAlleleCoverage cov2 = new HlaAlleleCoverage(allele2, 10,20, 0);
        HlaAlleleCoverage cov3 = new HlaAlleleCoverage(allele3, 10,20, 0);
        HlaAlleleCoverage cov4 = new HlaAlleleCoverage(allele4, 10,20, 0);
        HlaAlleleCoverage cov5 = new HlaAlleleCoverage(allele5, 10,20, 0);
        HlaAlleleCoverage cov6 = new HlaAlleleCoverage(allele6, 10,20, 0);
        List<HlaAlleleCoverage> coverageList = Lists.newArrayList(cov1, cov2, cov3, cov4, cov5, cov6);

        HlaComplexCoverage complexCoverage = HlaComplexCoverage.create(coverageList);
        assertEquals(60, complexCoverage.UniqueCoverage);
        assertEquals(180, complexCoverage.TotalCoverage);

        // homozygous solution
        coverageList = Lists.newArrayList(cov1, cov3, cov5);

        complexCoverage = HlaComplexCoverage.create(coverageList);
        complexCoverage.expandToSixAlleles();
        assertEquals(EXPECTED_ALLELE_COUNT, complexCoverage.getAlleleCoverage().size());
        assertEquals(30, complexCoverage.UniqueCoverage);
        assertEquals(90, complexCoverage.TotalCoverage);

        for(int i = 0; i < EXPECTED_ALLELE_COUNT; ++i)
        {
            assertEquals(15.0, complexCoverage.getAlleleCoverage().get(i).TotalCoverage, 0.01);
        }

        // test filling in missing coverage (applicable for tumor and RNA)
        coverageList = Lists.newArrayList(cov1, cov2, cov5);

        List<HlaAllele> alleleList = Lists.newArrayList(allele1, allele2, allele3, allele4, allele5, allele6);

        complexCoverage = HlaComplexCoverage.create(coverageList);
        complexCoverage.populateMissingCoverage(alleleList);
        assertEquals(EXPECTED_ALLELE_COUNT, complexCoverage.getAlleleCoverage().size());
        assertEquals(90, complexCoverage.TotalCoverage);

        assertEquals(30.0, complexCoverage.getAlleleCoverage().get(0).TotalCoverage, 0.01);
        assertEquals(30.0, complexCoverage.getAlleleCoverage().get(1).TotalCoverage, 0.01);
        assertEquals(0, complexCoverage.getAlleleCoverage().get(2).TotalCoverage, 0.01);
        assertEquals(0, complexCoverage.getAlleleCoverage().get(3).TotalCoverage, 0.01);
        assertEquals(30.0, complexCoverage.getAlleleCoverage().get(4).TotalCoverage, 0.01);
        assertEquals(0, complexCoverage.getAlleleCoverage().get(5).TotalCoverage, 0.01);

    }
}
