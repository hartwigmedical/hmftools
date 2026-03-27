package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;

import org.junit.Test;

public class HlaAlleleFactoryTest
{
    @Test
    public void testRegEx()
    {
        var matcher = HlaAlleleFactory.matchHlaRegEx("A*01:01");

        assertEquals("A", matcher.group("geneSymbol"));
        assertEquals("01", matcher.group("alleleGroup"));
        assertEquals("01", matcher.group("hlaProtein"));
    }

    @Test
    public void testMatchHlaRegExClassII()
    {
        var matcher = HlaAlleleFactory.matchHlaRegEx("DPB1*01:01");

        assertEquals("DPB1", matcher.group("geneSymbol"));
        assertEquals("01", matcher.group("alleleGroup"));
        assertEquals("01", matcher.group("hlaProtein"));
    }

    @Test
    public void testHlaAlleleConversion()
    {
        LilacRecord lilac = ImmutableLilacRecord.builder()
                .qc("PASS")
                .addAlleles(lilacAllele("A*02:01", 1.0, 100, 200, 50))
                .build();

        List<HlaAllele> result = HlaAlleleFactory.convertHlaAlleles(lilac, true, true);

        assertEquals(1, result.size());
        HlaAllele allele = result.get(0);
        assertEquals("HLA-A", allele.gene());
        assertEquals("02", allele.alleleGroup());
        assertEquals("01", allele.hlaProtein());
        assertEquals(HlaAllele.GeneClass.MHC_CLASS_I, allele.geneClass());
        assertEquals(1, allele.germlineCopyNumber());
        assertEquals(1.0, allele.tumorCopyNumber(), 0.01);
        assertEquals(Integer.valueOf(100), allele.refFragments());
        assertEquals(200, allele.tumorFragments());
        assertEquals(Integer.valueOf(50), allele.rnaFragments());
        assertEquals(Set.of(HlaAllele.QcStatus.PASS), allele.qcStatus());
        assertEquals("hlaAllele[A*02:01 1]", allele.findingKey());
        assertEquals("A*02:01", allele.allele());
    }

    @Test
    public void refAndRnaFragmentsAreNullWhenNotAvailable()
    {
        LilacRecord lilac = ImmutableLilacRecord.builder()
                .qc("PASS")
                .addAlleles(lilacAllele("C*01:02", 2.0, 50, 150, 30))
                .build();

        List<HlaAllele> result = HlaAlleleFactory.convertHlaAlleles(lilac, false, false);

        HlaAllele allele = result.get(0);
        assertNull(allele.refFragments());
        assertNull(allele.rnaFragments());
    }

    private static LilacAllele lilacAllele(String allele, double tumorCopyNumber, Integer refFragments, int tumorFragments,
            Integer rnaFragments)
    {
        return ImmutableLilacAllele.builder()
                .allele(allele)
                .tumorCopyNumber(tumorCopyNumber)
                .refFragments(refFragments)
                .tumorFragments(tumorFragments)
                .rnaFragments(rnaFragments)
                .somaticMissense(0.0)
                .somaticNonsenseOrFrameshift(0.0)
                .somaticSplice(0.0)
                .somaticSynonymous(0.0)
                .somaticInframeIndel(0.0)
                .build();
    }
}
