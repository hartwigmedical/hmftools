package com.hartwig.hmftools.viridian.profiling.extract;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;

public class ViralKmerIndexTest
{
    private static final String REFERENCE = "ACGATCGGATTCCAAGGTTCAACGTGACCTTAGGCATCGGATACCGTTAC";
    private static final String REFERENCE_WINDOW = REFERENCE.substring(10, 40);

    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testReadSharingReferenceSubsequenceIsDetected() throws IOException
    {
        assertTrue(index(REFERENCE).hasViralKmer(REFERENCE_WINDOW.getBytes()));
    }

    @Test
    public void testUnrelatedReadIsRejected() throws IOException
    {
        assertFalse(index(REFERENCE).hasViralKmer("TTTTGGGGCCCCAAAATTTTGGGG".getBytes()));
    }

    @Test
    public void testLowercaseReadMatchesUppercaseReference() throws IOException
    {
        assertTrue(index(REFERENCE).hasViralKmer(REFERENCE_WINDOW.toLowerCase().getBytes()));
    }

    @Test
    public void testReverseComplementReadIsDetected() throws IOException
    {
        String reverseComplement = Nucleotides.reverseComplementBases(REFERENCE_WINDOW);
        assertTrue(index(REFERENCE).hasViralKmer(reverseComplement.getBytes()));
    }

    @Test
    public void testReadShorterThanKmerLengthIsRejected() throws IOException
    {
        assertFalse(index(REFERENCE).hasViralKmer("ACGATCGGATTCCAAG".getBytes()));
    }

    // The reference's only k-mer spans an N, so nothing is indexed: a read with a T in place of the N is not detected.
    @Test
    public void testKmerSpanningAmbiguousReferenceBaseIsNotIndexed() throws IOException
    {
        assertFalse(index("ACGATCGGANTTCCAAGGT").hasViralKmer("ACGATCGGATTTCCAAGGT".getBytes()));
    }

    private ViralKmerIndex index(String... sequences) throws IOException
    {
        Path fasta = mTempDir.getRoot().toPath().resolve("virus_reference.fasta");
        List<String> fastaLines = new ArrayList<>();
        List<SAMSequenceRecord> dictionaryRecords = new ArrayList<>();
        for(int i = 0; i < sequences.length; ++i)
        {
            fastaLines.add(">v" + i);
            fastaLines.add(sequences[i]);
            dictionaryRecords.add(new SAMSequenceRecord("v" + i, sequences[i].length()));
        }
        Files.write(fasta, fastaLines);
        FastaSequenceIndexCreator.create(fasta, false);

        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(dictionaryRecords));
        StringWriter dictionary = new StringWriter();
        new SAMTextHeaderCodec().encode(dictionary, header);
        Files.writeString(fasta.resolveSibling("virus_reference.dict"), dictionary.toString());

        return ViralKmerIndex.build(fasta.toString());
    }
}
