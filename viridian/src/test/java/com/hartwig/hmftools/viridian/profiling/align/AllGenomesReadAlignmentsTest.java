package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.readAlignment;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.virusReference;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;

public class AllGenomesReadAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, new Virus("Virus v1", new VirusGroup("Group A")));

    @Test
    public void testFromOriginClipped()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(divergentAlignment("r1", 0), originClipped("r2"), originClipped("r3")));

        assertEquals(Set.of(ReadId.parse("r1")), allAlignments.store().byRead().keySet());
        assertEquals(Map.of(V1, 2), allAlignments.originClippedReads());
    }

    @Test
    public void testFromOriginClippedReadCount()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(originClipped("r1"), originClipped("r1"), originClipped("r2")));

        assertEquals(Map.of(V1, 2), allAlignments.originClippedReads());
    }

    @Test
    public void testFromRepeatAlignments()
    {
        // The aligner reports every placement it finds, so one read can be aligned to a virus genome several times. Only
        // its best alignment is kept, wherever it comes in the input, and how many it had is measured instead.

        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(
                        divergentAlignment("r1", 9),
                        divergentAlignment("r1", 4),
                        divergentAlignment("r1", 7),
                        divergentAlignment("r2", 1)));

        Map<ReadId, ViralReadAlignment> expected = Map.of(
                ReadId.parse("r1"), divergentAlignment("r1", 4),
                ReadId.parse("r2"), divergentAlignment("r2", 1));
        assertEquals(expected, allAlignments.store().byGenome().get(V1));
        assertEquals(Map.of(V1, List.of(3, 1)), allAlignments.alignmentCountsByGenome());
    }

    @Test
    public void testDecodeDivergence()
    {
        // Divergence is NM 1 plus 40 clipped bases. The MD tag puts the mismatch 30 bases after the start.

        ViralReadAlignment expected = new ViralReadAlignment(
                ReadId.parse("r1"), V1, 30, TextCigarCodec.decode("40S60M"), 60, 41, new int[] { 60 });

        assertEquals(expected, decode(record(30, "40S60M", 60, 1, "30A29")));
    }

    @Test
    public void testDecodeMissingTag()
    {
        assertThrows(IllegalStateException.class, () -> decode(record(30, "60M", null, 0, "60")));
        assertThrows(IllegalStateException.class, () -> decode(record(30, "60M", 60, null, "60")));
        assertThrows(IllegalStateException.class, () -> decode(record(30, "60M", 60, 0, null)));
    }

    private static ViralReadAlignment divergentAlignment(String readName, int divergence)
    {
        return readAlignment(readName, V1, 1, LENGTH + "M", LENGTH - divergence, divergence);
    }

    // A right clip projecting well past the virus genome end.
    private static ViralReadAlignment originClipped(String readName)
    {
        return readAlignment(readName, V1, 1, LENGTH + "M50S", LENGTH, 50);
    }

    private static ViralReadAlignment decode(SAMRecord record)
    {
        return new AllGenomesReadAlignments.AlignmentDecoder(virusReference(V1)).decode(record);
    }

    private static SAMRecord record(
            int start, String cigar, @Nullable Integer alignerScore, @Nullable Integer editDistance,
            @Nullable String mismatchesAndDeletions)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(V1.contigName(), V1.length()));

        SAMRecord record = new SAMRecord(header);
        record.setReadName("r1");
        record.setReferenceName(V1.contigName());
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases("A".repeat(TextCigarCodec.decode(cigar).getReadLength()).getBytes());
        if(alignerScore != null)
        {
            record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignerScore.intValue());
        }
        if(editDistance != null)
        {
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, editDistance.intValue());
        }
        if(mismatchesAndDeletions != null)
        {
            record.setAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE, mismatchesAndDeletions);
        }
        return record;
    }
}
