package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SaTagRewriterTest
{
    @Test
    public void testSupplementaryAdmissionRequiresALiftableSaAlignment()
    {
        ContigTranslator translator = new ContigTranslator(List.of(threeExonContig()));

        assertTrue(SaTagRewriter.hasLiftableAlignment(TX_CONTIG + ",1,+,50M,60,2;", translator));
        assertFalse(SaTagRewriter.hasLiftableAlignment(TX_CONTIG + ",9999,+,50M,60,2;", translator));
        assertFalse(SaTagRewriter.hasLiftableAlignment(null, translator));
    }

    @Test
    public void testPrimarySaUsesFinalEmittedSupplementaryFields()
    {
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M");
        SAMRecord supp = supplementaryRecord(CHR_1, 250, "30M20H", "stale");
        supp.setReadNegativeStrandFlag(true);
        supp.setMappingQuality(42);
        supp.setAttribute("NM", 3);

        String saTag = SaTagRewriter.buildForRecord(
                List.of(primary, supp), new boolean[] { true, true }, primary, 0);

        assertEquals(CHR_1 + ",250,-,30M20S,42,3;", saTag);
    }

    @Test
    public void testSupplementarySaUsesFinalPrimaryThenEmittedSiblings()
    {
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M");
        primary.setMappingQuality(55);
        primary.setAttribute("NM", 1);
        SAMRecord target = supplementaryRecord(CHR_1, 250, "30M20S", "stale");
        SAMRecord sibling = supplementaryRecord(CHR_1, 400, "20H30M", "stale");
        sibling.setReadNegativeStrandFlag(true);
        sibling.setMappingQuality(41);
        sibling.setAttribute("NM", 4);
        SAMRecord dropped = supplementaryRecord(CHR_1, 800, "25M25S", "stale");
        List<SAMRecord> records = List.of(primary, target, sibling, dropped);

        String saTag = SaTagRewriter.buildForRecord(
                records, new boolean[] { true, true, true, false }, primary, 1);

        assertEquals(
                CHR_1 + ",100,+,50M,55,1;" + CHR_1 + ",400,-,20S30M,41,4;",
                saTag);
    }

    @Test
    public void testPrimaryHasNoSaWhenNoSupplementaryWillBeEmitted()
    {
        SAMRecord primary = primaryRecord(CHR_1, 100, "50M");
        SAMRecord dropped = supplementaryRecord(CHR_1, 250, "30M20S", "stale");
        SAMRecord unmapped = supplementaryRecord(CHR_1, 400, "20S30M", "stale");
        unmapped.setReadUnmappedFlag(true);

        assertNull(SaTagRewriter.buildForRecord(
                List.of(primary, dropped, unmapped), new boolean[] { true, false, true }, primary, 0));
    }
}
