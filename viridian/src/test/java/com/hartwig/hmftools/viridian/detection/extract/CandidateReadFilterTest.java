package com.hartwig.hmftools.viridian.detection.extract;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.MateState.MATE_MAPPED_TO_HOST;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.MateState.MATE_MAPPED_TO_VIRUS;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.MateState.MATE_UNMAPPED;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.MateState.UNPAIRED;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.MAPPED_TO_HOST;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.MAPPED_TO_HOST_CLIPPED;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.MAPPED_TO_VIRUS;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.MAPPED_TO_VIRUS_CLIPPED;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.REDUX_UNMAPPED_FROM_HOST;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.REDUX_UNMAPPED_FROM_VIRUS;
import static com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilterTest.ReadState.REGULAR_UNMAPPED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class CandidateReadFilterTest
{
    private static final SAMFileHeader HEADER = new SAMFileHeader(new SAMSequenceDictionary(List.of(
            new SAMSequenceRecord("chr1", 10000),
            new SAMSequenceRecord("chrEBV", 171823))));

    private static final String BASES = "A".repeat(100);
    private static final int POSITION = 100;

    private static final CandidateReadFilter FILTER =
            new CandidateReadFilter(CandidateReadSource.ALL, 20, Set.of("chrEBV"), null);

    private static final CandidateReadFilter FULLY_UNMAPPED_AND_VIRUS_FILTER =
            new CandidateReadFilter(CandidateReadSource.FULLY_UNMAPPED_AND_VIRUS, 20, Set.of("chrEBV"), null);

    private static final boolean ACCEPT = true;
    private static final boolean REJECT = false;

    enum ReadState
    {
        MAPPED_TO_VIRUS("chrEBV", "100M", null),
        MAPPED_TO_VIRUS_CLIPPED("chrEBV", "30S70M", null),
        MAPPED_TO_HOST("chr1", "100M", null),
        MAPPED_TO_HOST_CLIPPED("chr1", "30S70M", null),
        REGULAR_UNMAPPED(null, "*", null),
        REDUX_UNMAPPED_FROM_VIRUS(null, "*", "chrEBV:500"),
        REDUX_UNMAPPED_FROM_HOST(null, "*", "chr1:500");

        @Nullable
        final String MappedContig;
        final String Cigar;
        @Nullable
        final String ReduxUnmapCoords;

        ReadState(@Nullable String mappedContig, String cigar, @Nullable String reduxUnmapCoords)
        {
            MappedContig = mappedContig;
            Cigar = cigar;
            ReduxUnmapCoords = reduxUnmapCoords;
        }
    }

    enum MateState
    {
        UNPAIRED(false, null),
        MATE_MAPPED_TO_VIRUS(true, "chrEBV"),
        MATE_MAPPED_TO_HOST(true, "chr1"),
        MATE_UNMAPPED(true, null);

        final boolean IsPaired;
        @Nullable
        final String MappedContig;

        MateState(boolean isPaired, @Nullable String mappedContig)
        {
            IsPaired = isPaired;
            MappedContig = mappedContig;
        }
    }

    record ReadMateExpectation(ReadState read, MateState mate, boolean isCandidateForAll, boolean isCandidateForFullyUnmappedAndVirus)
    {
    }

    // All plausible states of read and mate, with the expected filter result.
    // Easier to test programmatically than writing a test case for each.
    private static final List<ReadMateExpectation> READ_MATE_EXPECTATIONS = List.of(
            new ReadMateExpectation(MAPPED_TO_HOST, UNPAIRED, REJECT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST, MATE_MAPPED_TO_VIRUS, ACCEPT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST, MATE_MAPPED_TO_HOST, REJECT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST, MATE_UNMAPPED, ACCEPT, REJECT),

            new ReadMateExpectation(MAPPED_TO_HOST_CLIPPED, UNPAIRED, ACCEPT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST_CLIPPED, MATE_MAPPED_TO_VIRUS, ACCEPT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST_CLIPPED, MATE_MAPPED_TO_HOST, ACCEPT, REJECT),
            new ReadMateExpectation(MAPPED_TO_HOST_CLIPPED, MATE_UNMAPPED, ACCEPT, REJECT),

            new ReadMateExpectation(MAPPED_TO_VIRUS, UNPAIRED, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS, MATE_MAPPED_TO_VIRUS, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS, MATE_MAPPED_TO_HOST, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS, MATE_UNMAPPED, ACCEPT, ACCEPT),

            new ReadMateExpectation(MAPPED_TO_VIRUS_CLIPPED, UNPAIRED, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS_CLIPPED, MATE_MAPPED_TO_VIRUS, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS_CLIPPED, MATE_MAPPED_TO_HOST, ACCEPT, ACCEPT),
            new ReadMateExpectation(MAPPED_TO_VIRUS_CLIPPED, MATE_UNMAPPED, ACCEPT, ACCEPT),

            new ReadMateExpectation(REGULAR_UNMAPPED, UNPAIRED, ACCEPT, ACCEPT),
            new ReadMateExpectation(REGULAR_UNMAPPED, MATE_MAPPED_TO_VIRUS, ACCEPT, ACCEPT),
            new ReadMateExpectation(REGULAR_UNMAPPED, MATE_MAPPED_TO_HOST, ACCEPT, REJECT),
            new ReadMateExpectation(REGULAR_UNMAPPED, MATE_UNMAPPED, ACCEPT, ACCEPT),

            new ReadMateExpectation(REDUX_UNMAPPED_FROM_HOST, UNPAIRED, REJECT, REJECT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_HOST, MATE_MAPPED_TO_VIRUS, ACCEPT, ACCEPT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_HOST, MATE_MAPPED_TO_HOST, REJECT, REJECT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_HOST, MATE_UNMAPPED, REJECT, REJECT),

            new ReadMateExpectation(REDUX_UNMAPPED_FROM_VIRUS, UNPAIRED, ACCEPT, ACCEPT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_VIRUS, MATE_MAPPED_TO_VIRUS, ACCEPT, ACCEPT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_VIRUS, MATE_MAPPED_TO_HOST, ACCEPT, REJECT),
            new ReadMateExpectation(REDUX_UNMAPPED_FROM_VIRUS, MATE_UNMAPPED, ACCEPT, ACCEPT)
    );

    record ClipExpectation(String cigar, @Nullable String supplementaryAlignments, boolean isCandidate)
    {
    }

    @Test
    public void testIsCandidateRecordForEveryReadAndMateState()
    {
        List<ReadMateExpectation> wrongOutcomes = READ_MATE_EXPECTATIONS.stream()
                .map(expectation -> actualOutcome(expectation.read(), expectation.mate()))
                .filter(actual -> !READ_MATE_EXPECTATIONS.contains(actual))
                .toList();

        assertEquals(List.of(), wrongOutcomes);
    }

    @Test
    public void testReadMateExpectationsCoverEveryStateOnce()
    {
        Set<List<Enum<?>>> states = READ_MATE_EXPECTATIONS.stream()
                .map(expectation -> List.<Enum<?>>of(expectation.read(), expectation.mate()))
                .collect(Collectors.toSet());

        assertEquals(ReadState.values().length * MateState.values().length, READ_MATE_EXPECTATIONS.size());
        assertEquals(READ_MATE_EXPECTATIONS.size(), states.size());
    }

    private static final String HOST_SUPPLEMENTARY = "chr1,200,+,70M30S,60,0;";
    private static final String VIRUS_SUPPLEMENTARY = "chrEBV,200,+,70M30S,60,0;";

    // Soft clips on an unpaired host read, so the clip is the only possible evidence. The threshold is 20 bases.
    private static final List<ClipExpectation> CLIP_EXPECTATIONS = List.of(
            new ClipExpectation("100M", null, REJECT),

            new ClipExpectation("19S81M", null, REJECT),
            new ClipExpectation("20S80M", null, ACCEPT),
            new ClipExpectation("30S70M", null, ACCEPT),

            new ClipExpectation("81M19S", null, REJECT),
            new ClipExpectation("80M20S", null, ACCEPT),
            new ClipExpectation("70M30S", null, ACCEPT),

            // Each side is judged on its own, so two short clips do not add up.
            new ClipExpectation("19S62M19S", null, REJECT),
            new ClipExpectation("20S61M19S", null, ACCEPT),
            new ClipExpectation("19S61M20S", null, ACCEPT),

            // Clipped bases aligned elsewhere count only if every supplementary is on a virus decoy contig.
            new ClipExpectation("30S70M", HOST_SUPPLEMENTARY, REJECT),
            new ClipExpectation("30S70M", VIRUS_SUPPLEMENTARY, ACCEPT),
            new ClipExpectation("30S70M", VIRUS_SUPPLEMENTARY + HOST_SUPPLEMENTARY, REJECT),
            new ClipExpectation("70M30S", HOST_SUPPLEMENTARY, REJECT),
            new ClipExpectation("70M30S", VIRUS_SUPPLEMENTARY, ACCEPT),

            // A virus supplementary does not make up for a clip below the threshold.
            new ClipExpectation("19S81M", VIRUS_SUPPLEMENTARY, REJECT)
    );

    @Test
    public void testIsCandidateRecordForEverySoftClip()
    {
        List<ClipExpectation> wrongOutcomes = CLIP_EXPECTATIONS.stream()
                .filter(expectation -> isCandidateClippedHostRead(expectation) != expectation.isCandidate())
                .toList();

        assertEquals(List.of(), wrongOutcomes);
    }

    @Test
    public void testIsCandidateContigAllSourceAcceptsEveryContig()
    {
        assertTrue(FILTER.isCandidateContig("chr1"));
        assertTrue(FILTER.isCandidateContig("chrEBV"));
    }

    @Test
    public void testIsCandidateContigFullyUnmappedAndVirusSourceAcceptsOnlyVirus()
    {
        assertFalse(FULLY_UNMAPPED_AND_VIRUS_FILTER.isCandidateContig("chr1"));
        assertTrue(FULLY_UNMAPPED_AND_VIRUS_FILTER.isCandidateContig("chrEBV"));
    }

    @Test
    public void testIsCandidateRecordMalformedReduxUnmapAttributeRejected()
    {
        SAMRecord record = record(REGULAR_UNMAPPED, UNPAIRED);
        record.setAttribute(UNMAP_ATTRIBUTE, "chr1");

        assertThrows(RuntimeException.class, () -> FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordRejectsDuplicate()
    {
        SAMRecord record = record(MAPPED_TO_VIRUS, UNPAIRED);
        record.setDuplicateReadFlag(true);

        assertFalse(FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordRejectsSecondary()
    {
        SAMRecord record = record(MAPPED_TO_VIRUS, UNPAIRED);
        record.setSecondaryAlignment(true);

        assertFalse(FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordRejectsSupplementary()
    {
        SAMRecord record = record(MAPPED_TO_VIRUS, UNPAIRED);
        record.setSupplementaryAlignmentFlag(true);

        assertFalse(FILTER.isCandidateRecord(record));
    }

    private static ReadMateExpectation actualOutcome(ReadState read, MateState mate)
    {
        SAMRecord record = record(read, mate);
        return new ReadMateExpectation(
                read, mate, FILTER.isCandidateRecord(record), FULLY_UNMAPPED_AND_VIRUS_FILTER.isCandidateRecord(record));
    }

    private static boolean isCandidateClippedHostRead(ClipExpectation expectation)
    {
        SAMRecord record = record(MAPPED_TO_HOST, UNPAIRED);
        record.setCigarString(expectation.cigar());
        if(expectation.supplementaryAlignments() != null)
        {
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, expectation.supplementaryAlignments());
        }
        return FILTER.isCandidateRecord(record);
    }

    // Follows SAM placement: an unmapped read with a mapped mate takes its mate's position, and a mapped read records its
    // own position as that of its unmapped mate.
    private static SAMRecord record(ReadState read, MateState mate)
    {
        boolean isReadMapped = read.MappedContig != null;
        boolean isMateMapped = mate.MappedContig != null;

        SAMRecord record = new SAMRecord(HEADER);
        record.setReadName("read");
        record.setReadBases(BASES.getBytes());
        record.setCigarString(read.Cigar);

        record.setReadUnmappedFlag(!isReadMapped);
        String contig = isReadMapped ? read.MappedContig : isMateMapped ? mate.MappedContig : NO_ALIGNMENT_REFERENCE_NAME;
        setPosition(record, contig);

        record.setReadPairedFlag(mate.IsPaired);
        if(mate.IsPaired)
        {
            record.setFirstOfPairFlag(true);
            record.setMateUnmappedFlag(!isMateMapped);
            String mateContig = isMateMapped ? mate.MappedContig : isReadMapped ? read.MappedContig : NO_ALIGNMENT_REFERENCE_NAME;
            setMatePosition(record, mateContig);
        }

        if(read.ReduxUnmapCoords != null)
        {
            record.setAttribute(UNMAP_ATTRIBUTE, read.ReduxUnmapCoords);
        }
        return record;
    }

    private static void setPosition(SAMRecord record, String contig)
    {
        record.setReferenceName(contig);
        record.setAlignmentStart(contig.equals(NO_ALIGNMENT_REFERENCE_NAME) ? 0 : POSITION);
    }

    private static void setMatePosition(SAMRecord record, String contig)
    {
        record.setMateReferenceName(contig);
        record.setMateAlignmentStart(contig.equals(NO_ALIGNMENT_REFERENCE_NAME) ? 0 : POSITION);
    }
}
