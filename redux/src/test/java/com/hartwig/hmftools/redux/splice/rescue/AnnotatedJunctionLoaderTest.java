package com.hartwig.hmftools.redux.splice.rescue;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.util.Set;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class AnnotatedJunctionLoaderTest
{
    private File mDir;

    @Before
    public void setUp() throws IOException
    {
        mDir = Files.createTempDirectory("ensembl_loader_test").toFile();
    }

    @After
    public void tearDown()
    {
        if(mDir != null && mDir.exists())
        {
            for(File f : mDir.listFiles())
                f.delete();
            mDir.delete();
        }
    }

    private void writeFile(final String name, final String content) throws IOException
    {
        try(Writer w = new FileWriter(new File(mDir, name)))
        {
            w.write(content);
        }
    }

    @Test
    public void testLoadsIntronsFromAdjacentExonPairs() throws IOException
    {
        writeFile("ensembl_gene_data.csv",
                "GeneId,Chromosome\n"
                + "ENSG_TEST,1\n");

        // two transcripts, each with 3 exons producing 2 junctions
        writeFile("ensembl_trans_exon_data.csv",
                "GeneId,TransId,ExonStart,ExonEnd,ExonRank\n"
                + "ENSG_TEST,TX_A,1000,1099,1\n"
                + "ENSG_TEST,TX_A,1500,1599,2\n"
                + "ENSG_TEST,TX_A,2000,2099,3\n"
                + "ENSG_TEST,TX_B,1000,1099,1\n"     // shares first exon with TX_A
                + "ENSG_TEST,TX_B,1700,1799,2\n");   // different middle exon

        final Set<ChrIntron> introns = AnnotatedJunctionLoader.load(mDir.getAbsolutePath());

        // TX_A produces introns 1100-1499 and 1600-1999
        // TX_B produces intron 1100-1699
        assertEquals(3, introns.size());
        assertTrue(introns.contains(new ChrIntron("chr1", 1100, 1499)));
        assertTrue(introns.contains(new ChrIntron("chr1", 1600, 1999)));
        assertTrue(introns.contains(new ChrIntron("chr1", 1100, 1699)));
    }

    @Test
    public void testChrPrefixNormalization() throws IOException
    {
        writeFile("ensembl_gene_data.csv",
                "GeneId,Chromosome\n"
                + "ENSG_PREFIXED,chr3\n"
                + "ENSG_BARE,3\n");

        writeFile("ensembl_trans_exon_data.csv",
                "GeneId,TransId,ExonStart,ExonEnd,ExonRank\n"
                + "ENSG_PREFIXED,TX_P,100,199,1\n"
                + "ENSG_PREFIXED,TX_P,500,599,2\n"
                + "ENSG_BARE,TX_B,1000,1099,1\n"
                + "ENSG_BARE,TX_B,1500,1599,2\n");

        final Set<ChrIntron> introns = AnnotatedJunctionLoader.load(mDir.getAbsolutePath());

        // both should normalize to "chr3"
        assertEquals(2, introns.size());
        assertTrue(introns.contains(new ChrIntron("chr3", 200, 499)));
        assertTrue(introns.contains(new ChrIntron("chr3", 1100, 1499)));
    }

    @Test
    public void testSingleExonTranscriptsProduceNoIntrons() throws IOException
    {
        writeFile("ensembl_gene_data.csv",
                "GeneId,Chromosome\n"
                + "ENSG_X,1\n");

        writeFile("ensembl_trans_exon_data.csv",
                "GeneId,TransId,ExonStart,ExonEnd,ExonRank\n"
                + "ENSG_X,TX_SINGLE,1000,2000,1\n");

        final Set<ChrIntron> introns = AnnotatedJunctionLoader.load(mDir.getAbsolutePath());

        assertEquals(0, introns.size());
    }

    @Test
    public void testDuplicateJunctionsAcrossTranscriptsDeduplicate() throws IOException
    {
        writeFile("ensembl_gene_data.csv",
                "GeneId,Chromosome\n"
                + "ENSG_DUP,2\n");

        writeFile("ensembl_trans_exon_data.csv",
                "GeneId,TransId,ExonStart,ExonEnd,ExonRank\n"
                + "ENSG_DUP,T1,100,199,1\n"
                + "ENSG_DUP,T1,300,399,2\n"
                + "ENSG_DUP,T2,100,199,1\n"
                + "ENSG_DUP,T2,300,399,2\n");

        final Set<ChrIntron> introns = AnnotatedJunctionLoader.load(mDir.getAbsolutePath());

        // both transcripts produce the same intron 200-299; Set dedupes
        assertEquals(1, introns.size());
        assertTrue(introns.contains(new ChrIntron("chr2", 200, 299)));
    }

    @Test
    public void testReverseStrandRankOrder() throws IOException
    {
        // negative-strand transcripts have exons in higher-genome-pos → lower-genome-pos rank order.
        // loader should normalize so the intron still spans the lower exon end + 1 to higher start - 1.
        writeFile("ensembl_gene_data.csv",
                "GeneId,Chromosome\n"
                + "ENSG_NEG,4\n");

        writeFile("ensembl_trans_exon_data.csv",
                "GeneId,TransId,ExonStart,ExonEnd,ExonRank\n"
                + "ENSG_NEG,T_NEG,2000,2099,1\n"      // first by transcription = higher genome pos
                + "ENSG_NEG,T_NEG,1000,1099,2\n");    // second by transcription = lower genome pos

        final Set<ChrIntron> introns = AnnotatedJunctionLoader.load(mDir.getAbsolutePath());

        // intron spans the genomic gap regardless of rank ordering
        assertEquals(1, introns.size());
        assertTrue(introns.contains(new ChrIntron("chr4", 1100, 1999)));
    }
}
