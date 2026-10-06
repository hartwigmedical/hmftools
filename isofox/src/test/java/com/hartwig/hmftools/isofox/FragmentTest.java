package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.OTHER_TRANS;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNKNOWN;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNSPLICED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.isofox.common.Fragment;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentTest
{
    @Test
    public void testSingleEndFragment()
    {
        Read read = new Read(createRecord(100, "5S10M100N10M5S"));
        Fragment fragment = new Fragment(read);

        assertFalse(read.isReadPaired());
        assertEquals("fragment", fragment.id());
        assertEquals("1", fragment.chromosome());
        assertEquals(List.of(read), fragment.reads());
        assertThrows(UnsupportedOperationException.class, () -> fragment.reads().clear());
        assertEquals(1, fragment.fragmentCount());
        assertEquals(100, fragment.minAlignmentStart());
        assertEquals(219, fragment.maxAlignmentEnd());
        assertTrue(fragment.containsSplit());
        assertTrue(fragment.isFullyIntronic());
        assertTrue(fragment.uniqueValidRegions().isEmpty());
        assertFalse(fragment.spansMultipleRegions(1));
        assertFalse(fragment.readsInDifferentExons(1));
        assertEquals(List.of(new BaseRegion(100, 109), new BaseRegion(210, 219)), fragment.mergedMappings());
    }

    @Test
    public void testPairedFragmentInEitherOrder()
    {
        Read first = new Read(createRecord(200, "20M"));
        Read second = new Read(createRecord(100, "10S150M5S"));
        markAsPair(first, second);

        assertEquals(List.of(first, second), new Fragment(first, second).reads());
        assertEquals(List.of(second, first), new Fragment(second, first).reads());

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertEquals("fragment", fragment.id());
            assertEquals("1", fragment.chromosome());
            assertThrows(UnsupportedOperationException.class, () -> fragment.reads().add(first));
            assertEquals(1, fragment.fragmentCount());
            assertEquals(100, fragment.minAlignmentStart());
            assertEquals(249, fragment.maxAlignmentEnd());
            assertFalse(fragment.containsSplit());
            assertTrue(fragment.isFullyIntronic());
            assertEquals(List.of(new BaseRegion(100, 249)), fragment.mergedMappings());
        }
    }

    @Test
    public void testConsensusCountIsNotSummedAcrossMates()
    {
        SAMRecord firstRecord = createRecord(100, "20M");
        SAMRecord secondRecord = createRecord(200, "20M");
        firstRecord.setAttribute(CONSENSUS_READ_ATTRIBUTE, "5;5");
        secondRecord.setAttribute(CONSENSUS_READ_ATTRIBUTE, "5;5");
        Read first = new Read(firstRecord);
        Read second = new Read(secondRecord);

        assertEquals(5, new Fragment(first).fragmentCount());

        markAsPair(first, second);

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertEquals(5, fragment.fragmentCount());
        }
    }

    @Test
    public void testLeastAmbiguousReadDeterminesLoci()
    {
        Read unique = new Read(createRecord(100, "20M"));
        Read ambiguous = createMultiMappedRead("2,+500,20M,0;3,-700,20M,1;");
        Read lessAmbiguous = createMultiMappedRead("4,+900,20M,0;");

        assertEquals(1, new Fragment(unique).minNumLoci());
        assertNull(new Fragment(unique).altLoci());
        assertEquals(3, new Fragment(ambiguous).minNumLoci());
        assertEquals(2, new Fragment(ambiguous).altLoci().size());

        markAsPair(ambiguous, unique);

        for(Fragment fragment : List.of(new Fragment(ambiguous, unique), new Fragment(unique, ambiguous)))
        {
            assertEquals(1, fragment.minNumLoci());
            assertNull(fragment.altLoci());
        }

        markAsPair(ambiguous, lessAmbiguous);

        for(Fragment fragment : List.of(new Fragment(ambiguous, lessAmbiguous), new Fragment(lessAmbiguous, ambiguous)))
        {
            assertEquals(2, fragment.minNumLoci());
            assertEquals(1, fragment.altLoci().size());
            assertEquals(900, fragment.altLoci().get(0).Region.start());
        }
    }

    @Test
    public void testEqualLocusCountsKeepFirstSuppliedRead()
    {
        Read first = createMultiMappedRead("2,+500,20M,0;");
        Read second = createMultiMappedRead("3,+700,20M,0;");
        markAsPair(first, second);

        Fragment firstSupplied = new Fragment(first, second);
        Fragment secondSupplied = new Fragment(second, first);

        assertEquals(2, firstSupplied.minNumLoci());
        assertEquals(500, firstSupplied.altLoci().get(0).Region.start());
        assertEquals(2, secondSupplied.minNumLoci());
        assertEquals(700, secondSupplied.altLoci().get(0).Region.start());
    }

    @Test
    public void testMergedMappings()
    {
        Read split = new Read(createRecord(100, "10M100N10M"));
        Read overlapping = new Read(createRecord(105, "10M100N15M"));

        assertEquals(List.of(new BaseRegion(100, 109), new BaseRegion(210, 219)), new Fragment(split).mergedMappings());
        assertEquals(List.of(new BaseRegion(105, 114), new BaseRegion(215, 229)), new Fragment(overlapping).mergedMappings());

        markAsPair(split, overlapping);

        for(Fragment fragment : List.of(new Fragment(split, overlapping), new Fragment(overlapping, split)))
        {
            assertEquals(List.of(new BaseRegion(100, 114), new BaseRegion(210, 229)), fragment.mergedMappings());
        }

        Read splitRead = new Read(createRecord(100, "10M100N10M"));
        Read downstream = new Read(createRecord(300, "20M"));
        markAsPair(splitRead, downstream);

        for(Fragment fragment : List.of(new Fragment(splitRead, downstream), new Fragment(downstream, splitRead)))
        {
            assertTrue(fragment.containsSplit());
            assertEquals(List.of(new BaseRegion(100, 109), new BaseRegion(210, 219), new BaseRegion(300, 319)), fragment.mergedMappings());
        }

        Read upstream = new Read(createRecord(100, "30M"));
        Read adjacent = new Read(createRecord(120, "30M"));
        markAsPair(upstream, adjacent);

        for(Fragment fragment : List.of(new Fragment(upstream, adjacent), new Fragment(adjacent, upstream)))
        {
            assertEquals(List.of(new BaseRegion(100, 149)), fragment.mergedMappings());
        }
    }

    @Test
    public void testUniqueValidRegionsDeduplicatesAndSkipsExonIntron()
    {
        Read first = new Read(createRecord(100, "10M100N10M"));
        Read second = new Read(createRecord(100, "20M"));
        RegionReadData shared = createRegion("gene", 1, 1, "1", 100, 199);
        RegionReadData exonIntron = createRegion("gene", 1, 2, "1", 210, 299);
        RegionReadData unrelated = createRegion("other", 2, 1, "1", 300, 399);
        first.getMappedRegions().putAll(Map.of(shared, WITHIN_EXON, exonIntron, EXON_INTRON, unrelated, EXON_BOUNDARY));
        second.getMappedRegions().put(shared, EXON_BOUNDARY);

        assertFalse(new Fragment(first).isFullyIntronic());
        assertEquals(Set.of(shared, unrelated), Set.copyOf(new Fragment(first).uniqueValidRegions()));

        markAsPair(first, second);

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertFalse(fragment.isFullyIntronic());
            assertEquals(2, fragment.uniqueValidRegions().size());
            assertEquals(Set.of(shared, unrelated), Set.copyOf(fragment.uniqueValidRegions()));
        }
    }

    @Test
    public void testSpansMultipleRegionsCombinesBothReads()
    {
        RegionReadData exon1 = createRegion("gene", 1, 1, "1", 100, 199);
        RegionReadData exon2 = createRegion("gene", 1, 2, "1", 210, 299);

        Read exonIntronRead = new Read(createRecord(100, "10M100N10M"));
        exonIntronRead.getMappedRegions().putAll(Map.of(exon1, WITHIN_EXON, exon2, EXON_INTRON));
        assertFalse(new Fragment(exonIntronRead).spansMultipleRegions(1));

        Read bothExonsRead = new Read(createRecord(100, "10M100N10M"));
        bothExonsRead.getMappedRegions().putAll(Map.of(exon1, WITHIN_EXON, exon2, EXON_BOUNDARY));
        assertTrue(new Fragment(bothExonsRead).spansMultipleRegions(1));
        assertFalse(new Fragment(bothExonsRead).spansMultipleRegions(99));

        Read sameExonFirst = new Read(createRecord(100, "20M"));
        Read sameExonSecond = new Read(createRecord(120, "20M"));
        sameExonFirst.getMappedRegions().put(exon1, WITHIN_EXON);
        sameExonSecond.getMappedRegions().put(exon1, WITHIN_EXON);
        markAsPair(sameExonFirst, sameExonSecond);

        for(Fragment fragment : List.of(new Fragment(sameExonFirst, sameExonSecond), new Fragment(sameExonSecond, sameExonFirst)))
        {
            assertFalse(fragment.spansMultipleRegions(1));
        }

        Read exon1Read = new Read(createRecord(100, "20M"));
        Read exon2Read = new Read(createRecord(210, "20M"));
        exon1Read.getMappedRegions().put(exon1, WITHIN_EXON);
        exon2Read.getMappedRegions().put(exon2, EXON_BOUNDARY);
        assertFalse(new Fragment(exon1Read).spansMultipleRegions(1));

        markAsPair(exon1Read, exon2Read);

        for(Fragment fragment : List.of(new Fragment(exon1Read, exon2Read), new Fragment(exon2Read, exon1Read)))
        {
            assertTrue(fragment.spansMultipleRegions(1));
        }
    }

    @Test
    public void testPairedTranscriptSupportIsAnIntersection()
    {
        Read first = new Read(createRecord(100, "20M"));
        Read second = new Read(createRecord(200, "20M"));
        first.getTranscriptClassifications().putAll(Map.of(1, EXONIC, 2, SPLICE_JUNCTION, 3, ALT, 4, EXONIC, 5, EXONIC));
        second.getTranscriptClassifications().putAll(Map.of(1, SPLICE_JUNCTION, 2, EXONIC, 3, EXONIC, 5, ALT, 6, SPLICE_JUNCTION));

        assertEquals(Set.of(1, 2, 4, 5), Set.copyOf(new Fragment(first).validTypeTranscripts()));
        assertEquals(Set.of(3), new Fragment(first).invalidTranscripts(List.of(1, 2, 4, 5)));

        markAsPair(first, second);

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertEquals(2, fragment.validTypeTranscripts().size());
            assertEquals(Set.of(1, 2), Set.copyOf(fragment.validTypeTranscripts()));
            assertEquals(Set.of(3, 4, 5, 6), fragment.invalidTranscripts(List.of(1, 2)));
            assertEquals(Set.of(2, 3, 4, 5, 6), fragment.invalidTranscripts(List.of(1)));
            assertTrue(fragment.hasTranscriptClassification(1, EXONIC));
            assertTrue(fragment.hasTranscriptClassification(1, SPLICE_JUNCTION));
            assertTrue(fragment.hasTranscriptClassification(6, SPLICE_JUNCTION));
            assertFalse(fragment.hasTranscriptClassification(1, ALT));
            assertFalse(fragment.hasTranscriptClassification(99, EXONIC));
        }

        second.getTranscriptClassifications().clear();
        assertTrue(new Fragment(second).validTypeTranscripts().isEmpty());
        assertTrue(new Fragment(second).invalidTranscripts(List.of()).isEmpty());

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertTrue(fragment.validTypeTranscripts().isEmpty());
            assertEquals(Set.of(1, 2, 3, 4, 5), fragment.invalidTranscripts(List.of()));
        }
    }

    @Test
    public void testSetOtherTranscriptsPreservesAcceptedAndInvalidTypes()
    {
        for(boolean paired : List.of(false, true))
        {
            Read first = new Read(createRecord(100, "20M"));
            Read second = new Read(createRecord(200, "20M"));

            if(paired)
                markAsPair(first, second);

            Fragment fragment = paired ? new Fragment(first, second) : new Fragment(first);

            for(Read read : fragment.reads())
            {
                read.getTranscriptClassifications().putAll(Map.of(
                        1, EXONIC, 2, SPLICE_JUNCTION, 3, ALT, 4, UNSPLICED, 5, UNKNOWN, 6, OTHER_TRANS, 7, EXONIC));
            }

            fragment.setOtherTranscripts(List.of(1));

            for(Read read : fragment.reads())
            {
                assertEquals(Map.of(1, EXONIC, 2, OTHER_TRANS, 3, ALT, 4, UNSPLICED, 5, UNKNOWN, 6, OTHER_TRANS, 7, OTHER_TRANS),
                        read.getTranscriptClassifications());
            }

            assertEquals(List.of(1), fragment.validTypeTranscripts());
        }
    }

    @Test
    public void testDifferentExonsCompareRanksWithinTheRequestedTranscript()
    {
        Read first = new Read(createRecord(100, "20M"));
        Read second = new Read(createRecord(300, "20M"));
        first.getMappedRegions().put(createRegion("gene", 1, 1, "1", 100, 199), WITHIN_EXON);
        first.getMappedRegions().put(createRegion("other", 2, 2, "1", 200, 299), WITHIN_EXON);
        second.getMappedRegions().put(createRegion("gene", 1, 2, "1", 300, 399), WITHIN_EXON);

        assertFalse(new Fragment(first).readsInDifferentExons(1));

        markAsPair(first, second);

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertTrue(fragment.readsInDifferentExons(1));
        }

        second.getMappedRegions().put(createRegion("gene", 1, 1, "1", 100, 199), WITHIN_EXON);

        for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
        {
            assertFalse(fragment.readsInDifferentExons(1));
        }
    }

    @Test
    public void testOrientationUsesFirstOfPairRatherThanPositionOrListOrder()
    {
        for(boolean reversed : List.of(false, true))
        {
            Read first = new Read(createRecord(200, "20M"));
            Read second = new Read(createRecord(100, "20M"));
            first.bamRecord().setReadNegativeStrandFlag(reversed);
            second.bamRecord().setReadNegativeStrandFlag(!reversed);

            assertEquals(reversed ? REVERSE : FORWARD, new Fragment(first).orientation());

            markAsPair(first, second);

            for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
            {
                assertEquals(reversed ? REVERSE : FORWARD, fragment.orientation());
            }

            second.bamRecord().setReadNegativeStrandFlag(reversed);

            for(Fragment fragment : List.of(new Fragment(first, second), new Fragment(second, first)))
            {
                assertNull(fragment.orientation());
            }
        }
    }

    @Test
    public void testAdapterTrimmingRequiresAMate()
    {
        for(int order : List.of(0, 1))
        {
            Read first = new Read(createRecord(105, "10S20M10S"));
            Read second = new Read(createRecord(100, "10S20M10S"));
            second.bamRecord().setReadNegativeStrandFlag(true);
            String firstBases = first.readBases();
            String secondBases = second.readBases();

            for(Read read : List.of(first, second))
            {
                new Fragment(read).trimAdapterBases();
                assertEquals("10S20M10S", read.cigarStr());
                assertEquals(40, read.baseLength());
            }

            markAsPair(first, second);

            Fragment fragment = order == 0 ? new Fragment(first, second) : new Fragment(second, first);
            fragment.trimAdapterBases();

            assertEquals("10S20M5S", first.cigarStr());
            assertEquals("5S20M10S", second.cigarStr());
            assertEquals(firstBases.substring(0, 35), first.readBases());
            assertEquals(secondBases.substring(5), second.readBases());
            assertEquals(129, first.unclippedEnd());
            assertEquals(95, second.unclippedStart());
            assertEquals(List.of(new BaseRegion(100, 124)), fragment.mergedMappings());
        }
    }

    private static SAMRecord createRecord(int start, final String cigar)
    {
        int readLength = cigarFromStr(cigar).getReadLength();
        String bases = "ACGT".repeat((readLength + 3) / 4).substring(0, readLength);
        SAMRecord record = createSamRecord("fragment", "1", start, bases, cigar, "*", 0, false, false, null);
        record.setFlags(0);
        record.setMateAlignmentStart(0);
        record.setInferredInsertSize(0);
        return record;
    }

    private static Read createMultiMappedRead(final String xaTag)
    {
        SAMRecord record = createRecord(100, "20M");
        record.setAttribute(XA_ATTRIBUTE, xaTag);
        return new Read(record);
    }

    private static void markAsPair(final Read first, final Read second)
    {
        for(Read read : List.of(first, second))
        {
            Read mate = read == first ? second : first;
            SAMRecord record = read.bamRecord();
            record.setReadPairedFlag(true);
            record.setFirstOfPairFlag(read == first);
            record.setSecondOfPairFlag(read == second);
            record.setMateReferenceName(mate.chromosome());
            record.setMateAlignmentStart(mate.alignmentStart());
            record.setMateNegativeStrandFlag(mate.isReadReversed());
        }
    }
}
