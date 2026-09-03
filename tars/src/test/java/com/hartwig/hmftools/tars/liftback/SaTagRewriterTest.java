package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.liftedRecordAt;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.secondMateRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.supplementaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SaTagRewriterTest
{
    private static ContigTranslator newTranslator()
    {
        return new ContigTranslator(List.of(threeExonContig()));
    }

    @Test
    public void testRewriteSaTag()
    {
        ContigTranslator translator = newTranslator();

        // null or empty input -> null
        assertNull(SaTagRewriter.rewriteSaTag(null, translator));
        assertNull(SaTagRewriter.rewriteSaTag("", translator));

        // ref-contig entry passes through unchanged
        String refEntry = CHR_1 + ",1000,+,50M,60,2;";
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry, translator));

        // tx-contig entry lifted to genomic coordinates
        assertEquals(CHR_1 + ",100,+,50M,60,2;", SaTagRewriter.rewriteSaTag(TX_CONTIG + ",1,+,50M,60,2;", translator));

        // a malformed entry drops the whole tag: parsing is delegated to hmf-common's SupplementaryReadData, which
        // rejects the tag rather than skipping bad entries. bwa never emits malformed SA in practice.
        assertNull(SaTagRewriter.rewriteSaTag("junk;" + refEntry, translator));

        // duplicate lifted entries deduped
        assertEquals(refEntry, SaTagRewriter.rewriteSaTag(refEntry + refEntry, translator));

        // all entries fail -> null
        assertNull(SaTagRewriter.rewriteSaTag(TX_CONTIG + ",10000,+,50M,60,0;malformed;", translator));
    }

    @Test
    public void testExcludeKeyRemovesDroppedSuppAcrossClipType()
    {
        ContigTranslator translator = newTranslator();

        // The primary's SA lists a dropped supp soft-clipped, and the exclude key is soft too. The entry must be
        // removed rather than left dangling.
        Set<String> exclude = Set.of(CHR_1 + ":1000:-:19M247S");
        assertNull(SaTagRewriter.rewriteSaTag(CHR_1 + ",1000,-,19M247S,0,0;", translator, exclude));

        // and when the SA entry itself is hard-clipped, the H is normalised to S before matching the same key
        assertNull(SaTagRewriter.rewriteSaTag(CHR_1 + ",1000,-,19M247H,0,0;", translator, exclude));

        // an entry at a different locus is untouched by the exclude set
        assertEquals(
                CHR_1 + ",2000,-,19M247S,0,0;",
                SaTagRewriter.rewriteSaTag(CHR_1 + ",2000,-,19M247S,0,0;", translator, exclude));
    }

    @Test
    public void testDroppedSuppKeyMatchesTheSaEntryDescribingIt()
    {
        // The dropped supp's own RECORD is hard-clipped while the primary's SA entry describing it is soft-clipped.
        // droppedSuppSaKeys and rewriteSaTag have to key both the same way, or the primary keeps an SA entry
        // pointing at a record that was never emitted.
        ContigTranslator translator = newTranslator();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord droppedSupp = supplementaryRecord(TX_CONTIG, 1, "19M31H", TX_CONTIG + ",1,+,19M31S,60,0;");
        boolean[] willEmit = { true, false };

        Set<String> exclude = SaTagRewriter.droppedSuppSaKeys(
                List.of(primary, droppedSupp), willEmit, primary, translator);

        assertNull(SaTagRewriter.rewriteSaTag(TX_CONTIG + ",1,+,19M31S,60,0;", translator, exclude));
    }

    @Test
    public void testDroppedSuppSaKeysKeyedOffSuppBwaPlacement()
    {
        // A supp that will NOT be emitted contributes its own bwa placement lifted to the genome, not a post-liftback
        // trimmed cigar. Supp at tx-contig pos 1 (exon 1) lifts to chr1:100.
        ContigTranslator translator = newTranslator();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord droppedSupp = supplementaryRecord(TX_CONTIG, 1, "50M", TX_CONTIG + ",1,+,50M,0,0;");
        boolean[] willEmit = { true, false };

        Set<String> keys = SaTagRewriter.droppedSuppSaKeys(
                List.of(primary, droppedSupp), willEmit, primary, translator);

        assertEquals(Set.of(CHR_1 + ":100:+:50M"), keys);
    }

    @Test
    public void testDroppedSuppSaKeysSkipsEmittedSupp()
    {
        ContigTranslator translator = newTranslator();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 1, "50M", TX_CONTIG + ",1,+,50M,0,0;");
        boolean[] willEmit = { true, true };

        assertTrue(SaTagRewriter.droppedSuppSaKeys(List.of(primary, supp), willEmit, primary, translator).isEmpty());
    }

    @Test
    public void testDroppedSuppSaKeysSkipsPrimaryAndNonSupplementaryRecords()
    {
        // Only supplementaries other than the primary contribute. willEmit=false throughout isolates the
        // primary / supplementary-flag gate from the willEmit gate.
        ContigTranslator translator = newTranslator();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord mate = secondMateRecord(TX_CONTIG, 1, "50M");
        boolean[] willEmit = { false, false };

        assertTrue(SaTagRewriter.droppedSuppSaKeys(List.of(primary, mate), willEmit, primary, translator).isEmpty());
    }

    @Test
    public void testDroppedSuppSaKeysSkipsUnliftableSupp()
    {
        // A dropped supp whose bwa placement is an out-of-range tx position fails to lift -> no key emitted.
        ContigTranslator translator = newTranslator();
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord supp = supplementaryRecord(TX_CONTIG, 9999, "30M", TX_CONTIG + ",1,+,50M,0,0;");
        boolean[] willEmit = { true, false };

        assertTrue(SaTagRewriter.droppedSuppSaKeys(List.of(primary, supp), willEmit, primary, translator).isEmpty());
    }

    @Test
    public void testBuildSwappedSuppSaLeadsWithPrimaryThenSiblingSupps()
    {
        // records: [primary, supp being rebuilt (index 1), sibling supp (2)]. The rebuilt SA is the emitted primary
        // entry followed by the sibling's lifted entry; the record at suppIndex and the primary record are not
        // appended. Each entry renders chrom,pos,strand,cigar,mapQuality,0 with NM fixed at 0.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord suppRebuilt = supplementaryRecord(TX_CONTIG, 1, "30M20S", TX_CONTIG + ",1,+,50M,0,0;");
        SAMRecord siblingSupp = supplementaryRecord(TX_CONTIG, 1, "20S30M", TX_CONTIG + ",1,+,50M,0,0;");
        List<SAMRecord> records = List.of(primary, suppRebuilt, siblingSupp);

        LiftedRecord primaryResult = liftedRecordAt(CHR_1, 100, "50M", false);
        LiftedRecord siblingResult = liftedRecordAt(CHR_1, 200, "30M20S", true);
        LiftedRecord[] resolved = { primaryResult, null, siblingResult };
        boolean[] willEmit = { true, true, true };

        String saTag = SaTagRewriter.buildSwappedSuppSa(records, resolved, willEmit, primaryResult, 1);

        assertEquals(CHR_1 + ",100,+,50M,60,0;" + CHR_1 + ",200,-,30M20S,60,0;", saTag);
    }

    @Test
    public void testBuildSwappedSuppSaOmitsNonEmittedAndUnmappedSiblingSupps()
    {
        // Two sibling supps that must be excluded: one not emitting, one emitting but resolved to NO_ALIGNMENT_CIGAR.
        // The result is the primary entry alone.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord suppRebuilt = supplementaryRecord(TX_CONTIG, 1, "30M20S", TX_CONTIG + ",1,+,50M,0,0;");
        SAMRecord siblingNotEmitting = supplementaryRecord(TX_CONTIG, 1, "25M25S", TX_CONTIG + ",1,+,50M,0,0;");
        SAMRecord siblingUnmapped = supplementaryRecord(TX_CONTIG, 1, "20S30M", TX_CONTIG + ",1,+,50M,0,0;");
        List<SAMRecord> records = List.of(primary, suppRebuilt, siblingNotEmitting, siblingUnmapped);

        LiftedRecord primaryResult = liftedRecordAt(CHR_1, 100, "50M", false);
        LiftedRecord notEmittingResult = liftedRecordAt(CHR_1, 300, "25M25S", false);
        LiftedRecord unmappedResult = liftedRecordAt(CHR_1, 1, SAMRecord.NO_ALIGNMENT_CIGAR, false);
        LiftedRecord[] resolved = { primaryResult, null, notEmittingResult, unmappedResult };
        boolean[] willEmit = { true, true, false, true };

        String saTag = SaTagRewriter.buildSwappedSuppSa(records, resolved, willEmit, primaryResult, 1);

        assertEquals(CHR_1 + ",100,+,50M,60,0;", saTag);
    }

    @Test
    public void testBuildSwappedSuppSaNullWhenPrimaryUnplaced()
    {
        // Nothing lifted for the primary, so there is no anchor for the supp's SA.
        SAMRecord primary = primaryRecord(TX_CONTIG, 1, "50M");
        SAMRecord suppRebuilt = supplementaryRecord(TX_CONTIG, 1, "30M20S", TX_CONTIG + ",1,+,50M,0,0;");
        List<SAMRecord> records = List.of(primary, suppRebuilt);

        LiftedRecord primaryResult = LiftedRecord.unmapped("unliftable");
        LiftedRecord[] resolved = { primaryResult, null };
        boolean[] willEmit = { true, true };

        assertNull(SaTagRewriter.buildSwappedSuppSa(records, resolved, willEmit, primaryResult, 1));
    }
}
