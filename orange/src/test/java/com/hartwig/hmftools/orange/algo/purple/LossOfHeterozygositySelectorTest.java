package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LossOfHeterozygositySelectorTest
{
    private static final double EPSILON = 1.0E-10;

    @Test
    public void canSelectGenesForLOH()
    {
        String hrdGene = LossOfHeterozygositySelector.HRD_GENES.iterator().next();

        GeneCopyNumber hrdGeneWithLOH = GeneCopyNumberTestFactory.createGeneCopyNumber(hrdGene, 2D, 2D);
        GeneCopyNumber hrdGeneWithoutLOH = buildGeneCopyNumber(hrdGene, 2, 2, 2);

        String msiGene = LossOfHeterozygositySelector.MSI_GENES.iterator().next();
        GeneCopyNumber msiGeneWithLOH = buildGeneCopyNumber(msiGene, 2, 2, 0);
        GeneCopyNumber msiGeneWithoutLOH = buildGeneCopyNumber(msiGene, 0, 0, 0);

        GeneCopyNumber otherGeneWithLOH = buildGeneCopyNumber("other", 2, 2, 0);

        List<GeneCopyNumber> allGeneCopyNumbers = Lists.newArrayList(
                hrdGeneWithLOH, hrdGeneWithoutLOH, msiGeneWithLOH, msiGeneWithoutLOH, otherGeneWithLOH);

        List<GeneCopyNumber> all =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers, null, MicrosatelliteStatus.MSI, ChordStatus.HR_DEFICIENT);
        assertEquals(2, all.size());
        assertTrue(all.contains(hrdGeneWithLOH));
        assertTrue(all.contains(msiGeneWithLOH));

        List<GeneCopyNumber> msiOnly =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers, null, MicrosatelliteStatus.MSI, ChordStatus.HR_PROFICIENT);
        assertEquals(1, msiOnly.size());
        assertTrue(msiOnly.contains(msiGeneWithLOH));

        List<GeneCopyNumber> hrdOnly =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers, null, MicrosatelliteStatus.MSS, ChordStatus.HR_DEFICIENT);
        assertEquals(1, hrdOnly.size());
        assertTrue(hrdOnly.contains(hrdGeneWithLOH));

        List<GeneCopyNumber> none =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers, null, MicrosatelliteStatus.MSS, ChordStatus.HR_PROFICIENT);
        assertEquals(0, none.size());

        List<GeneCopyNumber> nullable =
                LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(allGeneCopyNumbers, null, MicrosatelliteStatus.MSS, null);
        assertEquals(0, nullable.size());
    }

    private static GeneCopyNumber buildGeneCopyNumber(
            final String gene, final double minCopyNumber, final double maxCopyNumber, final double minMinorAlleleCopyNumber)
    {
        return new GeneCopyNumber(CHR_1, 0, 0, gene, Strings.EMPTY, true,
                Strings.EMPTY, maxCopyNumber , minCopyNumber, minMinorAlleleCopyNumber, 1, 1,
                0, 0, 0, 1.0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN);
    }

    @Test
    public void canSelectGeneForLOHBasedOnGermlineDeletion()
    {
        String gene = LossOfHeterozygositySelector.HRD_GENES.iterator().next();
        GeneCopyNumber hrdGene = buildGeneCopyNumber(gene, 2, 2, 1);

        GermlineDeletion hetDeletion = GermlineDeletionTestFactory.create(gene, true, GermlineStatus.HET_DELETION, 1);
        List<GeneCopyNumber> lohList = runWithHRDOneGeneOneGermlineDeletion(hrdGene, hetDeletion);
        assertEquals(1, lohList.size());

        GeneCopyNumber loh = lohList.get(0);
        assertEquals(0, loh.MinMinorAlleleCopyNumber, EPSILON);
        assertEquals(1, loh.minCopyNumber(), EPSILON);

        GermlineDeletion homDeletion = GermlineDeletionTestFactory.create(gene, true, GermlineStatus.HOM_DELETION);
        assertTrue(runWithHRDOneGeneOneGermlineDeletion(hrdGene, homDeletion).isEmpty());
    }

    @NotNull
    private static List<GeneCopyNumber> runWithHRDOneGeneOneGermlineDeletion(@NotNull GeneCopyNumber geneCopyNumber,
            @NotNull GermlineDeletion germlineDeletion)
    {
        return LossOfHeterozygositySelector.selectHRDOrMSIGenesWithLOH(
                Lists.newArrayList(geneCopyNumber),
                Lists.newArrayList(germlineDeletion),
                MicrosatelliteStatus.MSS,
                ChordStatus.HR_DEFICIENT
        );
    }
}