package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class TransvalTest extends TransvalTestBase
{
    @Test
    public void mtorSNV()
    {
        // This example is based on TransvarConverterTest in the serve codebase.
        // The idea is to check that the current code agrees with Transval on
        // the key fields used by serve.
        TransvalSNV record = transval.calculateSNV("MTOR:p.L2230V");

        assertEquals("ENST00000361445", record.TranscriptId);
        assertEquals("1", record.Chromosome);
        assertEquals(11_122_101, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        assertEquals("A", record.ReferenceNucleotide);
        assertEquals("C", record.AlternateNucleotide);
        assertEquals("TTA", record.ReferenceCodon);
        assertEquals("GTA", record.AlternateCodons.get(0));
        assertEquals("GTC", record.AlternateCodons.get(1));
        assertEquals("GTG", record.AlternateCodons.get(2));
        assertEquals("GTT", record.AlternateCodons.get(3));
        assertEquals(4, record.AlternateCodons.size());
    }

    @Test
    public void brafSNV()
    {
    /*
    This test is based on the following log lines from a serve run:
21:35:18 - [DEBUG] - Running 'transvar panno --reference /data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --refversion hg38 --noheader --ensembl -i BRAF:p.V600E'
21:35:20 - [DEBUG] - Converting transvar output line to TransvarRecord: 'BRAF:p.V600E   ENST00000288602 (protein_coding)        BRAF    -       chr7:g.140753336A>T/c.1799T>A/p.V600E   inside_[cds_in_exon_15] CSQN=Missense;reference_codon=GTG;candidate_codons=GAG,GAA;candidate_mnv_variants=chr7:g.140753335_140753336delCAinsTT;aliases=ENSP00000288602;source=Ensembl'
21:35:20 - [DEBUG] - Interpreting transvar record: 'TransvarRecord{transcript=ENST00000288602, chromosome=7, gdnaPosition=140753336, variantSpanMultipleExons=false, annotation=TransvarSnvMnv{gdnaRef=A, gdnaAlt=T, referenceCodon=GTG, candidateCodons=[GAG, GAA]}}'
21:35:20 - [DEBUG] - Converted 'BRAF|null|p.V600E' to 2 hotspot(s)
21:35:20 - [INFO ] - Printing hotspots for 'BRAF:p.V600E' on transcript null
21:35:20 - [INFO ] -  Hotspot{ref=A, alt=T, chromosome=chr7, position=140753336}
21:35:20 - [INFO ] -  Hotspot{ref=CA, alt=TT, chromosome=chr7, position=140753335}
*/
        TransvalSNV record = transval.calculateSNV("BRAF:p.V600E");
//        assertEquals("ENST00000288602", record.TranscriptId); // Our ensembl data has ENST00000646891 as the canonical transcript
        assertEquals("7", record.Chromosome);
        assertEquals(140_753_336, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        assertEquals("A", record.ReferenceNucleotide);
        assertEquals("T", record.AlternateNucleotide);
        assertEquals("GTG", record.ReferenceCodon);
        assertEquals("GAG", record.AlternateCodons.get(0));
        assertEquals("GAA", record.AlternateCodons.get(1));
        assertEquals(2, record.AlternateCodons.size());
    }

    @Test
    public void vhl()
    {
        // Another example based on a test in serve's TransvarConverterTest.
        TransvalSNV record = transval.calculateSNV("VHL:p.G114R");
                assertEquals("ENST00000256474", record.TranscriptId);
        assertEquals("3", record.Chromosome);
//        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think
        assertTrue(record.SpansMultipleExons);
        assertEquals("G", record.ReferenceNucleotide);
        assertEquals("C", record.AlternateNucleotide);
        assertEquals("GGT", record.ReferenceCodon);
        // The serve test has the alternate codons in the order defined by Transvar,
        // which is simply the order in which they appear in its reverse codon table:
        // 'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
        // However, we order them by edit distance and then alphabetically.
        assertEquals("CGT", record.AlternateCodons.get(0));
        assertEquals("AGA", record.AlternateCodons.get(1));
        assertEquals("AGG", record.AlternateCodons.get(2));
        assertEquals("CGA", record.AlternateCodons.get(3));
        assertEquals("CGC", record.AlternateCodons.get(4));
        assertEquals("CGG", record.AlternateCodons.get(5));
        assertEquals(6, record.AlternateCodons.size());
    }
}
