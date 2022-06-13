package com.hartwig.hmftools.orange.algo.selection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.junit.Test;

public class LossOfHeterozygositySelectorTest {

    @Test
    public void canSelectGenesForLOH() {
        String hrdGene = LossOfHeterozygositySelector.HRD_GENES.iterator().next();
        GeneCopyNumber hrdGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName(hrdGene).minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();
        GeneCopyNumber hrdGeneWithoutLOH =
                GeneCopyNumberTestFactory.builder().geneName(hrdGene).minMinorAlleleCopyNumber(2D).minCopyNumber(2D).build();

        String msiGene = LossOfHeterozygositySelector.MSI_GENES.iterator().next();
        GeneCopyNumber msiGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName(msiGene).minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();
        GeneCopyNumber msiGeneWithoutLOH =
                GeneCopyNumberTestFactory.builder().geneName(msiGene).minMinorAlleleCopyNumber(0D).minCopyNumber(0D).build();

        GeneCopyNumber otherGeneWithLOH =
                GeneCopyNumberTestFactory.builder().geneName("other").minMinorAlleleCopyNumber(0D).minCopyNumber(2D).build();

        List<GeneCopyNumber> allGeneCopyNumbers =
                Lists.newArrayList(hrdGeneWithLOH, hrdGeneWithoutLOH, msiGeneWithLOH, msiGeneWithoutLOH, otherGeneWithLOH);

        List<GeneCopyNumber> all = LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSI,
                ChordStatus.HR_DEFICIENT);
        assertEquals(2, all.size());
        assertTrue(all.contains(hrdGeneWithLOH));
        assertTrue(all.contains(msiGeneWithLOH));

        List<GeneCopyNumber> msiOnly = LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSI,
                ChordStatus.HR_PROFICIENT);
        assertEquals(1, msiOnly.size());
        assertTrue(msiOnly.contains(msiGeneWithLOH));

        List<GeneCopyNumber> hrdOnly = LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSS,
                ChordStatus.HR_DEFICIENT);
        assertEquals(1, hrdOnly.size());
        assertTrue(hrdOnly.contains(hrdGeneWithLOH));

        List<GeneCopyNumber> none = LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers,
                MicrosatelliteStatus.MSS,
                ChordStatus.HR_PROFICIENT);
        assertEquals(0, none.size());
    }
}