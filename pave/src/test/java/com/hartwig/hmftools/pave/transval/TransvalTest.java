package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.junit.Test;

public class TransvalTest extends TransvalTestBase
{
    @Test
    public void mtorSNV()
    {
        // This example is based on TransvarConverterTest in the serve codebase.
        // The idea is to check that the current code agrees with Transval on
        // the key fields used by serve.
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("MTOR:p.L2230V");

        assertEquals("ENST00000361445", record.TranscriptId);
        assertEquals("1", record.Chromosome);
        assertEquals(11_122_101, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        assertEquals("A", record.ReferenceNucleotides);
        assertEquals("C", record.AlternateNucleotides);
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
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("BRAF:p.V600E");
//        assertEquals("ENST00000288602", record.TranscriptId); // Our ensembl data has ENST00000646891 as the canonical transcript
        assertEquals("7", record.Chromosome);
        assertEquals(140_753_336, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        assertEquals("A", record.ReferenceNucleotides);
        assertEquals("T", record.AlternateNucleotides);
        assertEquals("GTG", record.ReferenceCodon);
        assertEquals("GAG", record.AlternateCodons.get(0));
        assertEquals("GAA", record.AlternateCodons.get(1));
        assertEquals(2, record.AlternateCodons.size());
    }

    @Test
    public void vhlMultipleExons()
    {
        // Another example based on a test in serve's TransvarConverterTest.
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("VHL:p.G114R");
        assertEquals("ENST00000256474", record.TranscriptId);
        assertEquals("3", record.Chromosome);
//        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think
        assertTrue(record.SpansMultipleExons);
        assertEquals("G", record.ReferenceNucleotides);
        assertEquals("C", record.AlternateNucleotides);
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

    @Test
    public void tet2MNV()
    {
        // Another example based on a test in serve's TransvarConverterTest.
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("TET2:p.Y1294A");
        assertEquals("ENST00000380013", record.TranscriptId); // TransvarConvertTest has ENST00000540549
        assertEquals("4", record.Chromosome);
        //        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        assertEquals("TA", record.ReferenceNucleotides);
        assertEquals("GC", record.AlternateNucleotides);
        assertEquals("TAC", record.ReferenceCodon);

        assertEquals("GCC", record.AlternateCodons.get(0)); // Edit distance 2, the other options have edit distance 3
        assertEquals("GCA", record.AlternateCodons.get(1));
        assertEquals("GCG", record.AlternateCodons.get(2));
        assertEquals("GCT", record.AlternateCodons.get(3));
        assertEquals(4, record.AlternateCodons.size());
    }

    @Test
    public void egfrDelIns()
    {
/*
05:16:53 - [INFO ] - Printing hotspots for 'EGFR:p.747_750delinsP' on transcript null
05:16:53 - [INFO ] -  Hotspot{ref=TTAAGAGAAGCA, alt=CCT, chromosome=chr7, position=55174776}
05:16:53 - [INFO ] -  Hotspot{ref=TTAAGAGAAGCA, alt=CCG, chromosome=chr7, position=55174776}
05:16:53 - [INFO ] -  Hotspot{ref=TTAAGAGAAG, alt=C, chromosome=chr7, position=55174776}
05:16:53 - [INFO ] -  Hotspot{ref=TTAAGAGAAGCA, alt=CCC, chromosome=chr7, position=55174776}
 */
        TransvalInsertionDeletion record = (TransvalInsertionDeletion) transval.calculateVariant("EGFR:p.L747_A750delinsP");
        assertEquals("ENST00000275493", record.TranscriptId);
        assertEquals("7", record.Chromosome);
        assertEquals(55_174_776, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(10, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(4, hotspots.size());
        assertTrue(hotspots.contains(hotspot("TTAAGAGAAGCA", "CCT", "chr7", 55174776)));
        assertTrue(hotspots.contains(hotspot("TTAAGAGAAGCA", "CCG", "chr7", 55174776)));
        assertTrue(hotspots.contains(hotspot("TTAAGAGAAG", "C", "chr7", 55174776)));
        assertTrue(hotspots.contains(hotspot("TTAAGAGAAGCA", "CCC", "chr7", 55174776)));
    }

    @Test
    public void vhlDelIns()
    {
        /*
        Interpreting transvar record: 'TransvarRecord{
            transcript=ENST00000256474, chromosome=3, gdnaPosition=10141860, variantSpanMultipleExons=false,
            annotation=TransvarComplexInsertDelete{deletedBaseCount=12, insertedSequence=ATG, candidateAlternativeCodons=[ATG]}}'
        Converted 'VHL|null|p.A5_W8delinsM' to 1 hotspot(s)
        Printing hotspots for 'VHL:p.A5_W8delinsM' on transcript null
        Hotspot{ref=GCGGAGAACTG, alt=AT, chromosome=chr3, position=10141860}
         */
        TransvalInsertionDeletion record = (TransvalInsertionDeletion) transval.calculateVariant("VHL:p.A5_W8delinsM");
        assertEquals("ENST00000256474", record.TranscriptId);
        assertEquals("3", record.Chromosome);
        assertEquals(10_141_860, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(11, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(1, hotspots.size());
        assertTrue(hotspots.contains(hotspot("GCGGAGAACTG", "AT", "chr3", 10141860)));
    }

    @Test
    public void manyOptionsEGFR()
    {
        /*
        Converting transvar output line to TransvarRecord: 'EGFR:p.I744_K745delinsKIPVAI
        ENST00000275493 (protein_coding)	EGFR	+	chr7:g.55174767_55174772delinsAAGATCCCTGTAGCAATC/c.2230_2235delinsAAGATCCCTGTAGCAATC/p.I744_K745delinsKIPVAI
        inside_[cds_in_exon_19]	CSQN=MultiAAMissense;
        1152_CandidatesOmitted;aliases=ENSP00000275493;source=Ensembl'
Interpreting transvar record: 'TransvarRecord{transcript=ENST00000275493, chromosome=7, gdnaPosition=55174767, variantSpanMultipleExons=false, annotation=TransvarComplexInsertDelete{deletedBaseCount=6, insertedSequence=AAGATCCCTGTAGCAATC, candidateAlternativeCodons=[]}}'
Converted 'EGFR|null|p.I744_K745delinsKIPVAI' to 1 hotspot(s)
Printing hotspots for 'EGFR:p.I744_K745delinsKIPVAI' on transcript null
Hotspot{ref=TCAAG, alt=AGATCCCTGTAGCAATC, chromosome=chr7, position=55174768}
         */
        TransvalInsertionDeletion record = (TransvalInsertionDeletion) transval.calculateVariant("EGFR:p.I744_K745delinsKIPVAI");
        assertEquals("ENST00000275493", record.TranscriptId);
        assertEquals("7", record.Chromosome);
        assertEquals(55_174_768, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(5, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(1152, hotspots.size());
        assertTrue(hotspots.contains(hotspot("TCAAG", "AGATCCCTGTAGCAATC", "chr7", 55174768)));
    }

    @Test
    public void zyxDelInsAtEndOfExon()
    {
        /*
        Consider ZYX:p.P67_D70delinsWKY
        AAs 67-70 of ZYX are: P P E D
        Nucleotides are:     CCC CCG GAA G|AC
        The codon for D crosses from exon 1 into exon 2
        Options for W: {TGG}
        Options for K: {AAA, AAG}
        Options for Y, assuming that codon ends in AC: {TAC}
        So the coding change is: [CCC CCG GAA G] goes to [TGG AAA T] or [TGG AAG T]

        This is something that Transvar doesn't handle:
05:47:32 - [DEBUG] - Converted 'ZYX|null|p.P67_D70delinsWKY' to 0 hotspot(s)
05:47:32 - [INFO ] - Printing hotspots for 'ZYX:p.P67_D70delinsWKY' on transcript null
05:47:32 - [DEBUG] - Loaded 39500 genes from /data/resources/public/ensembl_data_cache/38/ensembl_gene_data.csv
05:47:34 - [DEBUG] - Loaded 39500 genes with 214977 transcripts and 1532157 exons from /data/resources/public/ensembl_data_cache/38/ensembl_trans_exon_data.csv
05:47:34 - [DEBUG] - Running 'transvar panno --reference /data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --refversion hg38 --noheader --ensembl -i ZYX:p.P67_D70delinsWKY'
05:47:36 - [DEBUG] - Converting transvar output line to TransvarRecord: 'ZYX:p.P67_D70delinsWKY	ENST00000322764 (protein_coding)	ZYX	+	chr7:g.143381770_143382249delinsTGGAAGTAT/c.199_210delinsTGGAAGTAT/p.P67_D70delinsWKY	inside_[cds_in_exons_[2,3]]	CSQN=MultiAAMissense;4_CandidatesOmitted;aliases=ENSP00000324422;source=Ensembl'
05:47:36 - [DEBUG] - Converting transvar output line to TransvarRecord: 'ZYX:p.P67_D70delinsWKY	ENST00000457235 (protein_coding)	ZYX	+	chr7:g.143381770_143382249delinsTGGAAGTAT/c.199_210delinsTGGAAGTAT/p.P67_D70delinsWKY	inside_[cds_in_exons_[1,2]]	CSQN=MultiAAMissense;4_CandidatesOmitted;aliases=ENSP00000400537;source=Ensembl'
05:47:36 - [DEBUG] - Interpreting transvar record: 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}}'
05:47:36 - [DEBUG] - Complex insert/delete spanning multiple exons. Ignoring 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}}'
05:47:36 - [WARN ] - Could not derive any hotspots from record TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}} for 'ZYX:p.P67_D70delinsWKY - null'
05:47:36 - [DEBUG] - Converted 'ZYX|null|p.P67_D70delinsWKY' to 0 hotspot(s)
05:47:36 - [INFO ] - Printing hotspots for 'ZYX:p.P67_D70delinsWKY' on transcript null
         */

        TransvalInsertionDeletion record = (TransvalInsertionDeletion) transval.calculateVariant("ZYX:p.P67_D70delinsWKY");
        assertEquals("ENST00000322764", record.TranscriptId);
        assertEquals("7", record.Chromosome);
        assertEquals(143_381_770, record.Position);
        assertTrue(record.SpansMultipleExons);

        assertEquals(10, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(hotspot("CCCCCGGAAG", "TGGAAAT", "chr7", 143_381_770)));
        assertTrue(hotspots.contains(hotspot("CCCCCGGAAG", "TGGAAGT", "chr7", 143_381_770)));
    }

    private TransvalHotspot hotspot(String ref, String alt, String chr, int position)
    {
        return new TransvalHotspot(ref, alt, chr, position);
    }
}
