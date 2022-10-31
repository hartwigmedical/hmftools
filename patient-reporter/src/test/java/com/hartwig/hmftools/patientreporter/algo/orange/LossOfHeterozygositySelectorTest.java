package com.hartwig.hmftools.patientreporter.algo.orange;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.patientreporter.algo.LohGenesReporting;
import com.hartwig.hmftools.patientreporter.util.Genes;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class LossOfHeterozygositySelectorTest {
    private static final double EPSILON = 1.0E-10;

    @Test
    public void canSelectGenesForLOH() {
        String hrdGene = Genes.HRD_GENES.iterator().next();
        GeneCopyNumber hrdGeneWithLOH = GeneCopyNumberTestFactory.builder()
                .geneName(hrdGene)
                .minMinorAlleleCopyNumber(0D)
                .minCopyNumber(2D)
                .chromosome("1")
                .chromosomeBand("p")
                .build();
        GeneCopyNumber hrdGeneWithoutLOH = GeneCopyNumberTestFactory.builder()
                .geneName(hrdGene)
                .minMinorAlleleCopyNumber(2D)
                .minCopyNumber(2D)
                .chromosome("1")
                .chromosomeBand("p")
                .build();

        String msiGene = Genes.MSI_GENES.iterator().next();
        GeneCopyNumber msiGeneWithLOH = GeneCopyNumberTestFactory.builder()
                .geneName(msiGene)
                .minMinorAlleleCopyNumber(0D)
                .minCopyNumber(2D)
                .chromosome("1")
                .chromosomeBand("p")
                .build();
        GeneCopyNumber msiGeneWithoutLOH = GeneCopyNumberTestFactory.builder()
                .geneName(msiGene)
                .minMinorAlleleCopyNumber(0D)
                .minCopyNumber(0D)
                .chromosome("1")
                .chromosomeBand("p")
                .build();

        GeneCopyNumber otherGeneWithLOH = GeneCopyNumberTestFactory.builder()
                .geneName("other")
                .minMinorAlleleCopyNumber(0D)
                .minCopyNumber(2D)
                .chromosome("1")
                .chromosomeBand("p")
                .build();

        List<GeneCopyNumber> allGeneCopyNumbers =
                Lists.newArrayList(hrdGeneWithLOH, hrdGeneWithoutLOH, msiGeneWithLOH, msiGeneWithoutLOH, otherGeneWithLOH);

        List<LohGenesReporting> msiOnly = LossOfHeterozygositySelector.selectMSIGenesWithLOH(allGeneCopyNumbers, MicrosatelliteStatus.MSI);
        assertEquals(1, msiOnly.size());
        assertEquals("MSH6", msiOnly.get(0).gene());
        assertEquals("1p", msiOnly.get(0).location());
        assertEquals(0D, msiOnly.get(0).minorAlleleCopies(), EPSILON);
        assertEquals(2D, msiOnly.get(0).tumorCopies(), EPSILON);

        List<LohGenesReporting> hrdOnly = LossOfHeterozygositySelector.selectHRDGenesWithLOH(allGeneCopyNumbers, ChordStatus.HR_DEFICIENT);
        assertEquals(1, hrdOnly.size());
        assertEquals("RAD51B", hrdOnly.get(0).gene());
        assertEquals("1p", hrdOnly.get(0).location());
        assertEquals(0D, msiOnly.get(0).minorAlleleCopies(), EPSILON);
        assertEquals(2D, msiOnly.get(0).tumorCopies(), EPSILON);

        List<LohGenesReporting> noneHRD = LossOfHeterozygositySelector.selectHRDGenesWithLOH(allGeneCopyNumbers, ChordStatus.HR_PROFICIENT);
        assertEquals(0, noneHRD.size());

        List<LohGenesReporting> noneMSI = LossOfHeterozygositySelector.selectMSIGenesWithLOH(allGeneCopyNumbers, MicrosatelliteStatus.MSS);
        assertEquals(0, noneMSI.size());
    }

}