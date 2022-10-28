package com.hartwig.hmftools.patientreporter.algo.orange;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.patientreporter.util.Genes;

import org.junit.Test;

public class LossOfHeterozygositySelectorTest {
    @Test
    public void canSelectGenesForLOH() {
        String hrdGene = Genes.HRD_GENES.iterator().next();
        GeneCopyNumber hrdGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName(hrdGene).minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();
        GeneCopyNumber hrdGeneWithoutLOH =
                GeneCopyNumberTestFactory.builder().geneName(hrdGene).minMinorAlleleCopyNumber(2D).minCopyNumber(2D).build();

        String msiGene = Genes.MSI_GENES.iterator().next();
        GeneCopyNumber msiGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName(msiGene).minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();
        GeneCopyNumber msiGeneWithoutLOH =
                GeneCopyNumberTestFactory.builder().geneName(msiGene).minMinorAlleleCopyNumber(0D).minCopyNumber(0D).build();

        GeneCopyNumber otherGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName("other").minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();

        List<GeneCopyNumber> allGeneCopyNumbers =
                Lists.newArrayList(hrdGeneWithLOH, hrdGeneWithoutLOH, msiGeneWithLOH, msiGeneWithoutLOH, otherGeneWithLOH);

        List<GeneCopyNumber> msiOnly = LossOfHeterozygositySelector.selectMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSI);
        assertEquals(1, msiOnly.size());
        assertTrue(msiOnly.contains(msiGeneWithLOH));

        List<GeneCopyNumber> hrdOnly = LossOfHeterozygositySelector.selectHRDGenesWithLOH(allGeneCopyNumbers,
                ChordStatus.HR_DEFICIENT);
        assertEquals(1, hrdOnly.size());
        assertTrue(hrdOnly.contains(hrdGeneWithLOH));

        List<GeneCopyNumber> noneHRD = LossOfHeterozygositySelector.selectHRDGenesWithLOH(allGeneCopyNumbers,
                ChordStatus.HR_PROFICIENT);
        assertEquals(0, noneHRD.size());

        List<GeneCopyNumber> noneMSI = LossOfHeterozygositySelector.selectMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSS);
        assertEquals(0, noneMSI.size());
    }

}