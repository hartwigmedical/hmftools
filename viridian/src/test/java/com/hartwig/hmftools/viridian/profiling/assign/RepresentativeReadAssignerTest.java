package com.hartwig.hmftools.viridian.profiling.assign;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.readAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignmentStore;
import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.jetbrains.annotations.Nullable;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RepresentativeReadAssignerTest
{
    @Rule
    public final TemporaryFolder mTempDir = new TemporaryFolder();

    private static final int READ_LENGTH = 100;

    private static final VirusGroup GROUP = new VirusGroup("Group A");
    private static final VirusGenome REPRESENTATIVE = new VirusGenome("v1", 1000, new Virus("Virus v1", GROUP));
    private static final VirusGenome TWIN = new VirusGenome("v2", 1000, new Virus("Virus v2", GROUP));
    private static final VirusGenome MINOR = new VirusGenome("v3", 1000, new Virus("Virus v3", GROUP));

    @Test
    public void testAssignReads()
    {
        ViralReadAlignmentStore alignments = new ViralReadAlignmentStore(List.of(
                alignment("r1/1", TWIN, 100, 3),
                alignment("r1/1", REPRESENTATIVE, 100, 7),
                alignment("r1/1", MINOR, 100, 1)));

        Map<ReadId, ViralReadAlignment> assignments =
                RepresentativeReadAssigner.assignReads(alignments.byRead(), Set.of(REPRESENTATIVE, TWIN));

        assertEquals(Map.of(ReadId.parse("r1/1"), alignment("r1/1", TWIN, 100, 3)), assignments);
    }

    @Test
    public void testAssignReadsNoRepresentativeAlignment()
    {
        // The read supports only a strain no group settled on, so it is evidence for nothing.

        ViralReadAlignmentStore alignments = new ViralReadAlignmentStore(List.of(alignment("r1/1", MINOR, 100, 1)));

        assertTrue(RepresentativeReadAssigner.assignReads(alignments.byRead(), Set.of(REPRESENTATIVE)).isEmpty());
    }

    @Test
    public void testAssignReadsTie()
    {
        // Equally good on both, so the choice must not depend on hit ordering.

        ViralReadAlignment onTwin = alignment("r1/1", TWIN, 100, 3);
        ViralReadAlignment onRepresentative = alignment("r1/1", REPRESENTATIVE, 100, 3);

        assertEquals(REPRESENTATIVE, assignedGenome(List.of(onTwin, onRepresentative)));
        assertEquals(REPRESENTATIVE, assignedGenome(List.of(onRepresentative, onTwin)));
    }

    @Test
    public void testAssignSharedStart() throws IOException
    {
        // One read with two alignments to the same genome from the same base, telling apart only by their CIGAR.
        // The chosen one must be copied, and only it.

        ViralReadAlignment chosen = readAlignment("r1/1", REPRESENTATIVE, 500, "53S38M21D49M8S", 40, 86);
        ViralReadAlignment other = readAlignment("r1/1", REPRESENTATIVE, 500, "95S53M", 42, 98);

        String sourceBam = writeBam(List.of(other, chosen), SAMFileHeader.SortOrder.coordinate);
        String outputBam = new File(mTempDir.getRoot(), "representative.bam").getPath();

        RepresentativeReadAssigner.assign(
                AllGenomesReadAlignments.from(Stream.of(other, chosen)).store(), Set.of(REPRESENTATIVE),
                sourceBam, outputBam);

        assertEquals(List.of("53S38M21D49M8S"), cigars(outputBam));
    }

    @Test
    public void testAssignOutputOrder() throws IOException
    {
        // The output is written presorted, relying on the coordinate sorted source. It must come out sorted and
        // indexed, with reads assigned across two representatives. r3 fits v1 better, so its v2 alignment is not copied.

        List<ViralReadAlignment> alignments = List.of(
                alignment("r1/1", TWIN, 300, 1),
                alignment("r2/1", REPRESENTATIVE, 50, 1),
                alignment("r3/1", REPRESENTATIVE, 700, 1),
                alignment("r3/1", TWIN, 10, 5),
                alignment("r4/1", TWIN, 20, 1));

        String sourceBam = writeBam(alignments, SAMFileHeader.SortOrder.coordinate);
        String outputBam = new File(mTempDir.getRoot(), "representative.bam").getPath();

        RepresentativeReadAssigner.assign(
                AllGenomesReadAlignments.from(alignments.stream()).store(), Set.of(REPRESENTATIVE, TWIN), sourceBam,
                outputBam);

        assertEquals(List.of("v1:50 r2/1", "v1:700 r3/1", "v2:20 r4/1", "v2:300 r1/1"), positions(outputBam, null));
        assertEquals(List.of("v2:20 r4/1", "v2:300 r1/1"), positions(outputBam, TWIN.contigName()));
    }

    @Test
    public void testAssignUnsortedSource() throws IOException
    {
        // The source BAM comes from the read aligner, which sorts. An unsorted one must fail, not write a broken BAM.

        List<ViralReadAlignment> alignments = List.of(
                alignment("r1/1", REPRESENTATIVE, 700, 1),
                alignment("r2/1", REPRESENTATIVE, 50, 1));

        String sourceBam = writeBam(alignments, SAMFileHeader.SortOrder.unsorted);
        String outputBam = new File(mTempDir.getRoot(), "representative.bam").getPath();
        ViralReadAlignmentStore store = AllGenomesReadAlignments.from(alignments.stream()).store();

        assertThrows(
                IllegalArgumentException.class,
                () -> RepresentativeReadAssigner.assign(store, Set.of(REPRESENTATIVE), sourceBam, outputBam));
    }

    private static VirusGenome assignedGenome(List<ViralReadAlignment> alignments)
    {
        ViralReadAlignmentStore store = new ViralReadAlignmentStore(alignments);
        Map<ReadId, ViralReadAlignment> assignments =
                RepresentativeReadAssigner.assignReads(store.byRead(), Set.of(REPRESENTATIVE, TWIN));
        return assignments.get(ReadId.parse("r1/1")).genome();
    }

    // The source BAM only needs the fields the copy is matched on, plus enough shape to be a valid record.
    // A coordinate sort order sorts the records. Unsorted keeps them in the given order.
    private String writeBam(List<ViralReadAlignment> alignments, SAMFileHeader.SortOrder sortOrder) throws IOException
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord(REPRESENTATIVE.contigName(), REPRESENTATIVE.length()),
                new SAMSequenceRecord(TWIN.contigName(), TWIN.length()),
                new SAMSequenceRecord(MINOR.contigName(), MINOR.length()))));
        header.setSortOrder(sortOrder);

        File bam = new File(mTempDir.getRoot(), "all.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, bam))
        {
            for(ViralReadAlignment alignment : alignments)
            {
                SAMRecord record = new SAMRecord(header);
                record.setReadName(alignment.readId().toString());
                record.setReferenceName(alignment.genome().contigName());
                record.setAlignmentStart(alignment.alignmentStart());
                record.setCigar(alignment.cigar());
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

    // "contig:start readName" per record, over the whole BAM or, through its index, one contig.
    private static List<String> positions(String bamFile, @Nullable String contig) throws IOException
    {
        List<String> positions = new ArrayList<>();
        try(SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFile));
                SAMRecordIterator records = contig == null ? reader.iterator() : reader.queryOverlapping(contig, 0, 0))
        {
            records.forEachRemaining(record -> positions.add(
                    record.getReferenceName() + ":" + record.getAlignmentStart() + " " + record.getReadName()));
        }
        return positions;
    }

    private static ViralReadAlignment alignment(String readName, VirusGenome genome, int start, int divergence)
    {
        return readAlignment(readName, genome, start, READ_LENGTH + "M", READ_LENGTH - divergence, divergence);
    }
}
