package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// covers LiftBackStats counter mechanics, composition / MAPQ tier derivation, and the summary TSV layout.
public class LiftBackStatsTest
{
    private Path mSummary;

    @Before
    public void setUp() throws IOException
    {
        mSummary = Files.createTempFile("liftback_summary_", ".tsv");
    }

    @After
    public void tearDown() throws IOException
    {
        Files.deleteIfExists(mSummary);
    }

    private static SAMRecord newRecord(final int mapq)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setMappingQuality(mapq);
        return record;
    }

    private static LiftBackResult resultWith(final LiftBackCategory category, final int numXaAlts,
            final List<LiftedAlignment> alignments)
    {
        return TarsTestFixtures.resultBuilder()
                .category(category)
                .comp(LiftBackResult.Composition.NONE)
                .pos(0).cigar("*")
                .inputMapq(0).updatedMapq(0)
                .numXaAlts(numXaAlts).numRefAlts(0)
                .numLoci(0).numDistinctCigarsAtPrimaryLocus(0)
                .refFullMatch(false)
                .alignments(alignments)
                .build();
    }

    private static LiftedAlignment refAlignment()
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "1", 100, "50M",
                "1", 100, "50M",
                0, 0,
                null, null, null,
                false, true);
    }

    private static LiftedAlignment txAlignment()
    {
        return new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "txContig", 1, "50M",
                "1", 100, "50M",
                0, 0,
                "ENST_X", "ENSG_X", "GENEX",
                false, true);
    }

    @Test
    public void testRecordsIncrementCountersByCategory()
    {
        final LiftBackStats stats = new LiftBackStats();
        stats.record(newRecord(60), resultWith(LiftBackCategory.REF_SINGLE, 0, List.of(refAlignment())));
        stats.record(newRecord(60), resultWith(LiftBackCategory.REF_SINGLE, 0, List.of(refAlignment())));
        stats.record(newRecord(0), resultWith(LiftBackCategory.UNMAPPED, 0, List.of()));

        assertEquals(3, stats.total());
        assertEquals(2, stats.categoryCount(LiftBackCategory.REF_SINGLE));
        assertEquals(1, stats.categoryCount(LiftBackCategory.UNMAPPED));
    }

    @Test
    public void testCompositionDerivation()
    {
        assertEquals(LiftBackResult.Composition.NONE,
                LiftBackResult.Composition.fromAlignments(List.of()));

        assertEquals(LiftBackResult.Composition.REF_ONLY,
                LiftBackResult.Composition.fromAlignments(List.of(refAlignment())));

        assertEquals(LiftBackResult.Composition.TX_ONLY,
                LiftBackResult.Composition.fromAlignments(List.of(txAlignment())));

        assertEquals(LiftBackResult.Composition.REF_AND_TX,
                LiftBackResult.Composition.fromAlignments(List.of(refAlignment(), txAlignment())));
    }

    @Test
    public void testMapqTierDerivation()
    {
        final LiftBackResult result = resultWith(LiftBackCategory.REF_SINGLE, 0, List.of(refAlignment()));

        assertEquals(LiftBackStats.MapqTier.MAPQ_ZERO,
                LiftBackStats.deriveMapqTier(newRecord(0), result));
        assertEquals(LiftBackStats.MapqTier.MAPQ_POS_UNIQUE,
                LiftBackStats.deriveMapqTier(newRecord(60), result));

        final LiftBackResult withAlts = resultWith(LiftBackCategory.REF_SINGLE, 1, List.of(refAlignment()));
        assertEquals(LiftBackStats.MapqTier.MAPQ_POS_MULTI,
                LiftBackStats.deriveMapqTier(newRecord(60), withAlts));
    }

    @Test
    public void testWriteSummaryEmitsExpectedTsvLayout() throws IOException
    {
        final LiftBackStats stats = new LiftBackStats();
        stats.record(newRecord(60), resultWith(LiftBackCategory.REF_SINGLE, 0, List.of(refAlignment())));
        stats.record(newRecord(0), resultWith(LiftBackCategory.UNMAPPED, 0, List.of()));

        stats.writeSummary(mSummary.toString());

        final List<String> lines = Files.readAllLines(mSummary);
        // header + at least one composition row + one category row
        assertTrue(lines.size() >= 3);
        assertEquals("Section\tRowKey\tColumnKey\tCount", lines.get(0));

        boolean hasRefOnlyMapqUnique = false;
        boolean hasARefGenome = false;
        for(final String line : lines)
        {
            if(line.equals("composition_x_mapq\tREF_ONLY\tMAPQ_POS_UNIQUE\t1"))
                hasRefOnlyMapqUnique = true;
            if(line.equals("category_x_mapq\tREF_SINGLE\tMAPQ_POS_UNIQUE\t1"))
                hasARefGenome = true;
        }
        assertTrue(hasRefOnlyMapqUnique);
        assertTrue(hasARefGenome);
    }
}
