package com.hartwig.hmftools.tars.liftback.supplementary;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

// Unit tests for AnnotatedJunctionIndex: exact-bounds contains() membership, the by-intron-start and by-intron-end
// lookups (keyed on chromosome+position, so several introns sharing one boundary all come back), and the empty /
// null-junction degenerate cases where every lookup must yield an empty (non-null) list.
public class AnnotatedJunctionIndexTest
{
    private static final String CHR_1 = "chr1";
    private static final String CHR_2 = "chr2";

    private static ChrBaseRegion intron(final String chromosome, final int start, final int end)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }

    private static AnnotatedJunctionIndex index(final ChrBaseRegion... introns)
    {
        return new AnnotatedJunctionIndex(new HashSet<>(List.of(introns)));
    }

    @Test
    public void testContainsMatchesIntronByValueNotIdentity()
    {
        // membership is value-based (ChrBaseRegion.equals), so a freshly built region with the same bounds hits.
        AnnotatedJunctionIndex junctionIndex = index(intron(CHR_1, 200, 299));

        assertTrue(junctionIndex.contains(intron(CHR_1, 200, 299)));
    }

    @Test
    public void testContainsRejectsRegionWithDifferentBounds()
    {
        // contains() is exact: a one-base end shift, or the same span on another chromosome, is not a member.
        AnnotatedJunctionIndex junctionIndex = index(intron(CHR_1, 200, 299));

        assertFalse(junctionIndex.contains(intron(CHR_1, 200, 300)));
        assertFalse(junctionIndex.contains(intron(CHR_2, 200, 299)));
    }

    @Test
    public void testIntroByStartReturnsAllIntronsSharingStart()
    {
        // two introns sharing a donor boundary (start 200) both come back keyed by that start.
        ChrBaseRegion shortIntron = intron(CHR_1, 200, 299);
        ChrBaseRegion longIntron = intron(CHR_1, 200, 450);
        AnnotatedJunctionIndex junctionIndex = index(shortIntron, longIntron);

        List<ChrBaseRegion> result = junctionIndex.introByStart(CHR_1, 200);

        assertEquals(2, result.size());
        assertTrue(result.contains(shortIntron));
        assertTrue(result.contains(longIntron));
    }

    @Test
    public void testIntroByEndReturnsAllIntronsSharingEnd()
    {
        // two introns sharing an acceptor boundary (end 299) both come back keyed by that end.
        ChrBaseRegion nearIntron = intron(CHR_1, 200, 299);
        ChrBaseRegion farIntron = intron(CHR_1, 50, 299);
        AnnotatedJunctionIndex junctionIndex = index(nearIntron, farIntron);

        List<ChrBaseRegion> result = junctionIndex.introByEnd(CHR_1, 299);

        assertEquals(2, result.size());
        assertTrue(result.contains(nearIntron));
        assertTrue(result.contains(farIntron));
    }

    @Test
    public void testIntroByStartMissReturnsEmptyList()
    {
        AnnotatedJunctionIndex junctionIndex = index(intron(CHR_1, 200, 299));

        assertTrue(junctionIndex.introByStart(CHR_1, 500).isEmpty());
    }

    @Test
    public void testIntroByEndMissReturnsEmptyList()
    {
        AnnotatedJunctionIndex junctionIndex = index(intron(CHR_1, 200, 299));

        assertTrue(junctionIndex.introByEnd(CHR_1, 500).isEmpty());
    }

    @Test
    public void testLookupsRequireMatchingChromosome()
    {
        // the boundary position matches but the chromosome does not, so neither lookup finds the intron.
        AnnotatedJunctionIndex junctionIndex = index(intron(CHR_1, 200, 299));

        assertTrue(junctionIndex.introByStart(CHR_2, 200).isEmpty());
        assertTrue(junctionIndex.introByEnd(CHR_2, 299).isEmpty());
    }

    @Test
    public void testEmptyIndexHasZeroSizeAndEmptyLookups()
    {
        AnnotatedJunctionIndex junctionIndex = index();

        assertEquals(0, junctionIndex.size());
        assertFalse(junctionIndex.contains(intron(CHR_1, 200, 299)));
        assertTrue(junctionIndex.introByStart(CHR_1, 200).isEmpty());
        assertTrue(junctionIndex.introByEnd(CHR_1, 299).isEmpty());
    }

    @Test
    public void testNullJunctionsTreatedAsEmptyIndex()
    {
        // the constructor null-guards the set, so a null argument behaves like an empty index rather than throwing.
        AnnotatedJunctionIndex junctionIndex = new AnnotatedJunctionIndex((Set<ChrBaseRegion>) null);

        assertEquals(0, junctionIndex.size());
        assertFalse(junctionIndex.contains(intron(CHR_1, 200, 299)));
        assertTrue(junctionIndex.introByStart(CHR_1, 200).isEmpty());
    }

    // ---- fromContigEntries: junctions derived from the sidecar's exon spans ----

    private static ContigEntry entry(final String chromosome, final int[]... exons)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] exon : exons)
        {
            spans.add(new BaseRegion(exon[0], exon[1]));
        }
        return ContigEntry.annotationOnly("g", "gn", "tn", chromosome, 1, spans);
    }

    @Test
    public void testIntronsFromAdjacentExonPairs()
    {
        // 3-exon entry -> 2 junctions; a second entry with a different middle exon -> 1 distinct junction
        AnnotatedJunctionIndex junctions = AnnotatedJunctionIndex.fromContigEntries(List.of(
                entry("chr1", new int[] { 1000, 1099 }, new int[] { 1500, 1599 }, new int[] { 2000, 2099 }),
                entry("chr1", new int[] { 1000, 1099 }, new int[] { 1700, 1799 })));

        assertEquals(3, junctions.size());
        assertTrue(junctions.contains(new ChrBaseRegion("chr1", 1100, 1499)));
        assertTrue(junctions.contains(new ChrBaseRegion("chr1", 1600, 1999)));
        assertTrue(junctions.contains(new ChrBaseRegion("chr1", 1100, 1699)));
    }

    @Test
    public void testSingleExonEntryProducesNoIntrons()
    {
        AnnotatedJunctionIndex junctions = AnnotatedJunctionIndex.fromContigEntries(List.of(
                entry("chr1", new int[] { 1000, 1999 })));

        assertEquals(0, junctions.size());
    }

    @Test
    public void testDuplicateJunctionsAcrossEntriesDeduplicate()
    {
        AnnotatedJunctionIndex junctions = AnnotatedJunctionIndex.fromContigEntries(List.of(
                entry("chr2", new int[] { 100, 199 }, new int[] { 300, 399 }),
                entry("chr2", new int[] { 100, 199 }, new int[] { 300, 399 })));

        assertEquals(1, junctions.size());
        assertTrue(junctions.contains(new ChrBaseRegion("chr2", 200, 299)));
    }

    @Test
    public void testExonsOutOfOrderStillGiveLowToHighIntron()
    {
        // spans are sorted by start, so exons supplied high-then-low still yield a low..high intron
        AnnotatedJunctionIndex junctions = AnnotatedJunctionIndex.fromContigEntries(List.of(
                entry("chr4", new int[] { 2000, 2099 }, new int[] { 1000, 1099 })));

        assertEquals(1, junctions.size());
        assertTrue(junctions.contains(new ChrBaseRegion("chr4", 1100, 1999)));
    }
}
