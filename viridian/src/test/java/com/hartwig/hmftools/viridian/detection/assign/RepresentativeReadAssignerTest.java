package com.hartwig.hmftools.viridian.detection.assign;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.detection.align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.align.AllAlignments;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RepresentativeReadAssignerTest
{
    @Rule
    public final TemporaryFolder mTempDir = new TemporaryFolder();

    private static final int READ_LENGTH = 100;

    private static final OncologyGroup GROUP = new OncologyGroup("Group A");
    private static final ViralContig REPRESENTATIVE = new ViralContig("v1", 1000, "Virus v1", GROUP);
    private static final ViralContig TWIN = new ViralContig("v2", 1000, "Virus v2", GROUP);
    private static final ViralContig MINOR = new ViralContig("v3", 1000, "Virus v3", GROUP);

    @Test
    public void testAssignReadsPicksLeastDivergentRepresentative()
    {
        ViralReadAlignments alignments = new ViralReadAlignments(List.of(
                alignment("r1/1", TWIN, 100, 3),
                alignment("r1/1", REPRESENTATIVE, 100, 7),
                alignment("r1/1", MINOR, 100, 1)));

        Map<ReadId, ViralReadAlignment> assignments =
                RepresentativeReadAssigner.assignReads(alignments.byRead(), Set.of(REPRESENTATIVE, TWIN));

        assertEquals(Map.of(ReadId.parse("r1/1"), alignment("r1/1", TWIN, 100, 3)), assignments);
    }

    // The read supports only a strain no group settled on, so it is evidence for nothing.
    @Test
    public void testAssignReadsDropsReadWithNoRepresentativeAlignment()
    {
        ViralReadAlignments alignments = new ViralReadAlignments(List.of(alignment("r1/1", MINOR, 100, 1)));

        assertTrue(RepresentativeReadAssigner.assignReads(alignments.byRead(), Set.of(REPRESENTATIVE)).isEmpty());
    }

    // Equally good on both, so the choice must not depend on hit ordering.
    @Test
    public void testAssignReadsBreaksTieByContig()
    {
        ViralReadAlignments alignments = new ViralReadAlignments(List.of(
                alignment("r1/1", TWIN, 100, 3),
                alignment("r1/1", REPRESENTATIVE, 100, 3)));

        Map<ReadId, ViralReadAlignment> assignments =
                RepresentativeReadAssigner.assignReads(alignments.byRead(), Set.of(REPRESENTATIVE, TWIN));

        assertEquals(REPRESENTATIVE, assignments.get(ReadId.parse("r1/1")).contig());
    }

    // One read with two alignments to the same contig from the same base, telling apart only by their CIGAR.
    // The chosen one must be copied, and only it.
    @Test
    public void testAssignCopiesTheChosenOfTwoAlignmentsSharingAStart() throws IOException
    {
        ViralReadAlignment chosen = new ViralReadAlignment(
                ReadId.parse("r1/1"), REPRESENTATIVE, 500, 586, "53S38M21D49M8S", 53, 8, 40, 86,
                List.of(new AlignedInterval(500, 38), new AlignedInterval(559, 49)));
        ViralReadAlignment other = new ViralReadAlignment(
                ReadId.parse("r1/1"), REPRESENTATIVE, 500, 552, "95S53M", 95, 0, 42, 98,
                List.of(new AlignedInterval(500, 53)));

        String sourceBam = writeBam("all.bam", List.of(other, chosen));
        String outputBam = new File(mTempDir.getRoot(), "representative.bam").getPath();

        RepresentativeReadAssigner.assign(
                AllAlignments.from(List.of(other, chosen), 148).alignments(), Set.of(REPRESENTATIVE),
                sourceBam, outputBam);

        assertEquals(List.of("53S38M21D49M8S"), cigars(outputBam));
    }

    // The source BAM only needs the fields the copy is matched on, plus enough shape to be a valid record.
    private String writeBam(String name, List<ViralReadAlignment> alignments) throws IOException
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord(REPRESENTATIVE.name(), REPRESENTATIVE.length()),
                new SAMSequenceRecord(TWIN.name(), TWIN.length()),
                new SAMSequenceRecord(MINOR.name(), MINOR.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        File bam = new File(mTempDir.getRoot(), name);
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, bam))
        {
            for(ViralReadAlignment alignment : alignments)
            {
                SAMRecord record = new SAMRecord(header);
                record.setReadName(alignment.readId().toString());
                record.setReferenceName(alignment.contig().name());
                record.setAlignmentStart(alignment.alignmentStart());
                record.setCigarString(alignment.cigar());
                record.setReadBases(SAMRecord.NULL_SEQUENCE);
                writer.addAlignment(record);
            }
        }
        return bam.getPath();
    }

    private static List<String> cigars(String bamFile) throws IOException
    {
        List<String> cigars = new ArrayList<>();
        try(SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile)))
        {
            reader.forEach(record -> cigars.add(record.getCigarString()));
        }
        return cigars;
    }

    private static ViralReadAlignment alignment(String readName, ViralContig contig, int start, int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, start, start + READ_LENGTH - 1, READ_LENGTH + "M", 0, 0, READ_LENGTH - divergence,
                divergence,
                List.of(new AlignedInterval(start, READ_LENGTH)));
    }
}
