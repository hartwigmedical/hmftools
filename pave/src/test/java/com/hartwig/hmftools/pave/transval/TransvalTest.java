package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.junit.Assert;
import org.junit.Test;

public class TransvalTest extends TransvalTestBase
{
    @Test
    public void mtorSNV()
    {
        // This example is based on TransvarConverterTest in the serve codebase.
        // The idea is to check that the current code agrees with Transval on
        // the key fields used by serve.
        /*
        05:47:12 - [DEBUG] - Converting transvar output line to TransvarRecord: 'MTOR:p.L2230V	ENST00000361445 (protein_coding)	MTOR	-	chr1:g.11122101A>C/c.6688T>G/p.L2230V	inside_[cds_in_exon_48]	CSQN=Missense;reference_codon=TTA;candidate_codons=GTA,GTC,GTG,GTT;candidate_mnv_variants=chr1:g.11122099_11122101delTAAinsGAC,chr1:g.11122099_11122101delTAAinsCAC,chr1:g.11122099_11122101delTAAinsAAC;aliases=ENSP00000354558;source=Ensembl'
05:47:12 - [DEBUG] - Interpreting transvar record: 'TransvarRecord{transcript=ENST00000361445, chromosome=1, gdnaPosition=11122101, variantSpanMultipleExons=false, annotation=TransvarSnvMnv{gdnaRef=A, gdnaAlt=C, referenceCodon=TTA, candidateCodons=[GTA, GTC, GTG, GTT]}}'
05:47:12 - [DEBUG] - Converted 'MTOR|null|p.L2230V' to 4 hotspot(s)
05:47:12 - [INFO ] - Printing hotspots for 'MTOR:p.L2230V' on transcript null
05:47:12 - [INFO ] -  Hotspot{ref=A, alt=C, chromosome=chr1, position=11122101}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=GAC, chromosome=chr1, position=11122099}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=CAC, chromosome=chr1, position=11122099}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=AAC, chromosome=chr1, position=11122099}
         */
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("MTOR:p.L2230V");

        assertEquals("ENST00000361445", record.transcriptId());
        assertEquals("1", record.Chromosome);
        //        assertEquals(11_122_101, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);
        //        assertEquals("A", record.ReferenceNucleotides);
        assertEquals("C", record.AlternateNucleotides);

        assertEquals(4, record.hotspots().size());
        assertTrue(record.hotspots().contains(hotspot("A", "C", "chr1", 11_122_101)));
        assertTrue(record.hotspots().contains(hotspot("TAA", "GAC", "chr1", 11_122_099)));
        assertTrue(record.hotspots().contains(hotspot("TAA", "CAC", "chr1", 11_122_099)));
        assertTrue(record.hotspots().contains(hotspot("TAA", "AAC", "chr1", 11_122_099)));

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
        //        assertEquals(140_753_336, record.Position); // serve example has 11_182_158, which is from v37, I think
        assertFalse(record.SpansMultipleExons);

        assertEquals(2, record.hotspots().size());
        assertTrue(record.hotspots().contains(hotspot("A", "T", "chr7", 140753336)));
        assertTrue(record.hotspots().contains(hotspot("CA", "TT", "chr7", 140753335)));

        //        assertEquals("A", record.ReferenceNucleotides);
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
        /*
        Converted 'VHL|null|p.G114R' to 1 hotspot(s)
        Printing hotspots for 'VHL:p.G114R' on transcript null
        Hotspot{ref=G, alt=C, chromosome=chr3, position=10142187}
         */
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("VHL:p.G114R");
        assertEquals("ENST00000256474", record.transcriptId());
        assertEquals("3", record.Chromosome);
        //        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think
        assertTrue(record.SpansMultipleExons);
        //        assertEquals("G", record.ReferenceNucleotides);
        assertEquals("C", record.AlternateNucleotides);
        assertEquals("GGT", record.ReferenceCodon);

        // todo - should we report all of the items below as hotspots?????????

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
        /*
        Converting transvar output line to TransvarRecord: 'TET2:p.Y1294A	ENST00000380013 (protein_coding)	TET2	+	chr4:g.105259695_105259696delTAinsGC/c.3880_3881delTAinsGC/p.Y1294A	inside_[cds_in_exon_7]	CSQN=Missense;reference_codon=TAC;candidate_codons=GCA,GCC,GCG,GCT;candidate_mnv_variants=chr4:g.105259695_105259697delTACinsGCA,chr4:g.105259695_105259697delTACinsGCG,chr4:g.105259695_105259697delTACinsGCT;aliases=ENSP00000369351;source=Ensembl'
        Interpreting transvar record: 'TransvarRecord{transcript=ENST00000380013, chromosome=4, gdnaPosition=105259695, variantSpanMultipleExons=false, annotation=TransvarSnvMnv{gdnaRef=TA, gdnaAlt=GC, referenceCodon=TAC, candidateCodons=[GCA, GCC, GCG, GCT]}}'
        Converted 'TET2|null|p.Y1294A' to 4 hotspot(s)
        Printing hotspots for 'TET2:p.Y1294A' on transcript null
        Hotspot{ref=TAC, alt=GCA, chromosome=chr4, position=105259695}
        Hotspot{ref=TA, alt=GC, chromosome=chr4, position=105259695}
        Hotspot{ref=TAC, alt=GCG, chromosome=chr4, position=105259695}
        Hotspot{ref=TAC, alt=GCT, chromosome=chr4, position=105259695}
         */
        TransvalSnvMnv record = (TransvalSnvMnv) transval.calculateVariant("TET2:p.Y1294A");
        assertEquals("ENST00000380013", record.transcriptId()); // TransvarConvertTest has ENST00000540549
        assertEquals("4", record.Chromosome);
        //        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think
        assertFalse(record.SpansMultipleExons);

        assertEquals(4, record.hotspots().size());
        assertTrue(record.hotspots().contains(hotspot("TA", "GC", "chr4", 105_259_695)));
        assertTrue(record.hotspots().contains(hotspot("TAC", "GCA", "chr4", 105_259_695)));
        assertTrue(record.hotspots().contains(hotspot("TAC", "GCG", "chr4", 105_259_695)));
        assertTrue(record.hotspots().contains(hotspot("TAC", "GCT", "chr4", 105_259_695)));

        //        assertEquals("TA", record.ReferenceNucleotides);
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
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("EGFR:p.L747_A750delinsP");
        assertEquals("ENST00000275493", record.transcriptId());
        assertEquals("7", record.Chromosome);
        //        assertEquals(55_174_776, record.Position);
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
    public void egfrServeExample()
    {
        /*
        Here's what Transvar says:
        Printing hotspots for 'EGFR:p.A750_I759delinsG' on transcript null
        Hotspot{ref=CAACATCTCCGAAAGCCAACAAGGAAATC, alt=GT, chromosome=chr7, position=55174786}
        Hotspot{ref=CAACATCTCCGAAAGCCAACAAGGAAATC, alt=GG, chromosome=chr7, position=55174786}
        Hotspot{ref=CAACATCTCCGAAAGCCAACAAGGAAATC, alt=GA, chromosome=chr7, position=55174786}
        Hotspot{ref=GCAACATCTCCGAAAGCCAACAAGGAAAT, alt=GG, chromosome=chr7, position=55174785}

        Here's my analysis:
        Section of transcript: ATSPKANKEI
        Corresponding nukes: GCAACATCTCCGAAAGCCAACAAGGAAATC
        A   T   S   P   K   A   N   K   E   I
        GCA ACA TCT CCG AAA GCC AAC AAG GAA ATC

        candidateAlternativeCodons={GGT, GGC, GGA, GGG}
        EXON starts at 55_174_722
        GGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAA - these 63 nukes precede the codon for A in the exon
        therefore A is at 55_174_722 + 63 = 55_174_785
        Changes keep the G of the codon for A, so begin at 55_174_785 + 1
         */
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("EGFR:p.A750_I759delinsG");
        assertEquals("ENST00000275493", record.transcriptId());
        assertEquals("7", record.Chromosome);
        //        assertEquals(55_174_786, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(28, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(4, hotspots.size());
        assertEquals(4, hotspots.size());
        assertTrue(hotspots.contains(hotspot("CAACATCTCCGAAAGCCAACAAGGAAATC", "GT", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(hotspot("CAACATCTCCGAAAGCCAACAAGGAAATC", "GG", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(hotspot("CAACATCTCCGAAAGCCAACAAGGAAATC", "GA", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(hotspot("CAACATCTCCGAAAGCCAACAAGGAAAT", "G", "chr7", 55_174_786)));
    }

    @Test
    public void egfrServeExample2()
    {
        /*
        This is an example where Transvar reports one hotspot, and we report 48.
        L747_K754delinsSPQ

        Here's what Transvar says:
        Printing hotspots for 'EGFR:p.L747_K754delinsSPQ' on transcript null
        Hotspot{ref=TTAAGAGAAGCAACATCTCCGA, alt=AGCCCTC, chromosome=chr7, position=55174776}

        L   R   E   A   T   S   P   K
        TTA AGA GAA GCA ACA TCT CCG AAA
        EXON starts at 55_174_722
        GGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAA - these 54 nukes precede the codon for L
        therefore L is at 55_174_722 + 54 = 55_174_776

        candidateAlternativeCodons={TCA, AGC, AGT, TCC, TCG, TCT}*{CCA, CCC, CCG, CCT}*{CAA, CAG} (48 possibilities)
         */
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("EGFR:p.L747_K754delinsSPQ");
        assertEquals("ENST00000275493", record.transcriptId());
        assertEquals("7", record.Chromosome);
        assertFalse(record.SpansMultipleExons);

        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(48, hotspots.size());
        assertTrue(hotspots.contains(hotspot("TTAAGAGAAGCAACATCTCCGA", "AGCCCTC", "chr7", 55_174_776)));
        String originalBases = StringUtils.remove("TTA AGA GAA GCA ACA TCT CCG AAA", ' ');
        for(TransvalHotspot hotspot : hotspots)
        {
            String newBases = StringUtils.replaceOnce(originalBases, hotspot.Ref, hotspot.Alt);
            AminoAcidSequence newAminoAcids = AminoAcidSequence.fromNucleotides(newBases);
            Assert.assertEquals("SPQ", newAminoAcids.toString());
        }
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
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("VHL:p.A5_W8delinsM");
        assertEquals("ENST00000256474", record.transcriptId());
        assertEquals("3", record.Chromosome);
        //        assertEquals(10_141_860, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(11, record.deletedBasesCount());
        checkSingleHotspot(record, "GCGGAGAACTG", "AT", "chr3", 10141860);
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
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("EGFR:p.I744_K745delinsKIPVAI");
        assertEquals("ENST00000275493", record.transcriptId());
        assertEquals("7", record.Chromosome);
        //        assertEquals(55_174_768, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(5, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(1152, hotspots.size());
        assertTrue(hotspots.contains(hotspot("TCAAG", "AGATCCCTGTAGCAATC", "chr7", 55174768)));
    }

    @Test
    public void deletionInsertionAtEndOfExon()
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
        Converted 'ZYX|null|p.P67_D70delinsWKY' to 0 hotspot(s)
        Printing hotspots for 'ZYX:p.P67_D70delinsWKY' on transcript null
        Loaded 39500 genes from /data/resources/public/ensembl_data_cache/38/ensembl_gene_data.csv
        Loaded 39500 genes with 214977 transcripts and 1532157 exons from /data/resources/public/ensembl_data_cache/38/ensembl_trans_exon_data.csv
        Running 'transvar panno --reference /data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --refversion hg38 --noheader --ensembl -i ZYX:p.P67_D70delinsWKY'
        Converting transvar output line to TransvarRecord: 'ZYX:p.P67_D70delinsWKY	ENST00000322764 (protein_coding)	ZYX	+	chr7:g.143381770_143382249delinsTGGAAGTAT/c.199_210delinsTGGAAGTAT/p.P67_D70delinsWKY	inside_[cds_in_exons_[2,3]]	CSQN=MultiAAMissense;4_CandidatesOmitted;aliases=ENSP00000324422;source=Ensembl'
        Converting transvar output line to TransvarRecord: 'ZYX:p.P67_D70delinsWKY	ENST00000457235 (protein_coding)	ZYX	+	chr7:g.143381770_143382249delinsTGGAAGTAT/c.199_210delinsTGGAAGTAT/p.P67_D70delinsWKY	inside_[cds_in_exons_[1,2]]	CSQN=MultiAAMissense;4_CandidatesOmitted;aliases=ENSP00000400537;source=Ensembl'
        Interpreting transvar record: 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}}'
        Complex insert/delete spanning multiple exons. Ignoring 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}}'
        Could not derive any hotspots from record TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143381770, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=480, insertedSequence=TGGAAGTAT, candidateAlternativeCodons=[]}} for 'ZYX:p.P67_D70delinsWKY - null'
        Converted 'ZYX|null|p.P67_D70delinsWKY' to 0 hotspot(s)
        Printing hotspots for 'ZYX:p.P67_D70delinsWKY' on transcript null
         */

        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("ZYX:p.P67_D70delinsWKY");
        assertEquals("ENST00000322764", record.transcriptId());
        assertEquals("7", record.Chromosome);
        //        assertEquals(143_381_770, record.Position);
        assertTrue(record.SpansMultipleExons);

        assertEquals(10, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(hotspot("CCCCCGGAAG", "TGGAAAT", "chr7", 143_381_770)));
        assertTrue(hotspots.contains(hotspot("CCCCCGGAAG", "TGGAAGT", "chr7", 143_381_770)));
    }

    @Test
    public void deletionInsertionAtStartOfExon()
    {
        // See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS5883
        // At position
        // 350 + LTLKEVEELEQLTQQLMQDMEHPQRQNVAVN.length
        // There is an E coded as G|AA split across exons 5 and 6.
        // The following AAs are L C G coded as CTC TGC GGC
        // Changing the start of exon 6 from
        // AA CTC TGC GGC
        // to
        // {AT, AC} {CCT, CCC, CCA, CCG}
        // would give
        // E L C G => D P
        // Note that exon 6 starts at 143_388_489
        //  G|AA CTC TGC GGC -> G|AC CCC: del A CTC TGC GG ins CCC

        /*
        03:00:58 - [DEBUG] - Converting transvar output line to TransvarRecord: 'ZYX:p.E382_G385delinsDP	ENST00000322764 (protein_coding)	ZYX	+	chr7:g.143388339_143388499delinsGACCCT/c.1144_1155delinsGACCCT/p.E382_G385delinsDP	inside_[cds_in_exons_[6,7]]	CSQN=MultiAAMissense;candidate_alternative_sequence=GAC/GAT+CCT/CCG/CCA/CCC;aliases=ENSP00000324422;source=Ensembl'
03:00:58 - [DEBUG] - Interpreting transvar record: 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143388339, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=161, insertedSequence=GACCCT, candidateAlternativeCodons=[GAC, GAT+CCT, CCG, CCA, CCC]}}'
03:00:58 - [DEBUG] - Complex insert/delete spanning multiple exons. Ignoring 'TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143388339, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=161, insertedSequence=GACCCT, candidateAlternativeCodons=[GAC, GAT+CCT, CCG, CCA, CCC]}}'
03:00:58 - [WARN ] - Could not derive any hotspots from record TransvarRecord{transcript=ENST00000322764, chromosome=7, gdnaPosition=143388339, variantSpanMultipleExons=true, annotation=TransvarComplexInsertDelete{deletedBaseCount=161, insertedSequence=GACCCT, candidateAlternativeCodons=[GAC, GAT+CCT, CCG, CCA, CCC]}} for 'ZYX:p.E382_G385delinsDP - null'
03:00:58 - [DEBUG] - Converted 'ZYX|null|p.E382_G385delinsDP' to 0 hotspot(s)
03:00:58 - [INFO ] - Printing hotspots for 'ZYX:p.E382_G385delinsDP' on transcript null
         */

        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("ZYX:p.E382_G385delinsDP");
        assertEquals("ENST00000322764", record.transcriptId());
        assertEquals("7", record.Chromosome);
        //        assertEquals(143_388_490, record.Position);
        assertTrue(record.SpansMultipleExons);

        assertEquals(9, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(8, hotspots.size());
        // G|AA CTC TGC GGC -> G|AT CCC
        assertTrue(hotspots.contains(hotspot("ACTCTGCGG", "TCC", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "TCCT", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "TCCA", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "TCCG", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGG", "CCC", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "CCCT", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "CCCA", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(hotspot("ACTCTGCGGC", "CCCG", "chr7", 143_388_490)));
    }

    @Test
    public void deletionInsertionPositiveStrand()
    {
        /*
        VHL:
 10141848
        1   2   3   4   5   6   ...
        M   P   R   R   A   E   ...
        ATG CCC CGG AGG GCG GAG ...

        R3_R4delinsQ
 10141854
        CGG AGG -> CAG|CAA. Options are GGAG -> A and GGAGG -> AA
        NOTE: Transvar gives essentially the same results, but has a common 'C' prefix
        for the ref and alt of the second hotspot:
        Printing hotspots for 'VHL:p.R3_R4delinsQ' on transcript null
        Hotspot{ref=GGAGG, alt=AA, chromosome=chr3, position=10141855}
        Hotspot{ref=CGGAG, alt=CA, chromosome=chr3, position=10141854}
        */
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("VHL:p.R3_R4delinsQ");
        assertEquals("ENST00000256474", record.transcriptId());
        assertEquals("3", record.Chromosome);
        //        assertEquals(10_141_855, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(4, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(hotspot("GGAG", "A", "chr3", 10_141_855)));
        assertTrue(hotspots.contains(hotspot("GGAGG", "AA", "chr3", 10_141_855)));
    }

    @Test
    public void deletionInsertionNegativeStrand()
    {
        /*
        BRAF:
140924703
        1   2   3   4   5   6   ...
        M   A   A   L   S   G   ...
        ATG GCG GCG CTG AGC GGT ...

        A3_L4delinsE
140924697
        GCG CTG -> GAA|GAG. Simplest option: GCGCTG -> GAG, so CGCT -> A. C is 140924696

        Reverse strand
        CAG CGC -> TTC|CTC. Options: CAG CG -> CT|TT, equal complexity. Position is 140924696 - 4

        NOTE: Transvar gives essentially these results, but leaves a leading 'C' on both the ref
        and alt of its first hotspot:
        Printing hotspots for 'BRAF:p.A3_L4delinsE' on transcript null
        Hotspot{ref=CAGCG, alt=CT, chromosome=chr7, position=140924692}
        Hotspot{ref=CAGCG, alt=TT, chromosome=chr7, position=140924692}
         */
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("BRAF:p.A3_L4delinsE");
        assertEquals("ENST00000646891", record.transcriptId()); // canonical
        assertEquals("7", record.Chromosome);
        //        assertEquals(140_924_693, record.Position);
        assertFalse(record.SpansMultipleExons);

        assertEquals(4, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(hotspot("AGCG", "T", "chr7", 140_924_693)));
        assertTrue(hotspots.contains(hotspot("CAGCG", "TT", "chr7", 140_924_692)));
    }

    @Test
    public void pik3r1delinsTest()
    {
        // A bug found when testing against the serve data.
        /*
        PIK3R1 E451_Y452delinsD
        '+' strand, chr5, exon 10
        E   Y
        GAA TAT
        D: {GAC, GAT}
        Exon starts at 68_293_709
        These bases precede E:
        GATCAAGTTGTCAAAGAAGATAATATTGAAGCTGTAGGGAAAAAATTACAT (51 bases)
        So the codon for E is at 68_293_709 + 51 = 68_293_760

        Need GAA TAT -> GAC or GAA TAT -> GAT
        ATAT -> C or AAT -> ""
         */
        TransvalDeletionInsertion record = (TransvalDeletionInsertion) transval.calculateVariant("PIK3R1:p.E451_Y452delinsD");
        assertEquals("ENST00000521381", record.transcriptId()); // canonical
        assertEquals("5", record.Chromosome);
        //        assertEquals(140_924_693, record.Position);
        assertFalse(record.SpansMultipleExons);

        //        assertEquals(4, record.deletedBasesCount());
        Set<TransvalHotspot> hotspots = record.hotspots();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(hotspot("AAT", "", "chr5", 68_293_761)));
        assertTrue(hotspots.contains(hotspot("ATAT", "C", "chr5", 68_293_762)));
    }

    @Test
    public void delSingleAminoAcidTest()
    {
        TransvalVariant record = transval.calculateVariant("PIK3R1:p.K459del");
        // EK: GAA AAA, G is at 68_293_781
        assertEquals("ENST00000521381", record.transcriptId()); // canonical
        assertEquals("5", record.Chromosome);
        assertFalse(record.SpansMultipleExons);
        checkSingleHotspot(record, "GAAA", "G", "chr5", 68_293_781);

        record = transval.calculateVariant("PIK3R1:p.D464del");
        //YDR: TAT GAT AGA, G of the D is at 68_293_799
        checkHotspots(record, hotspot("TATG", "T", "chr5", 68_293_796));
    }

    @Test
    public void delTwoAminoAcidsTest()
    {
        TransvalVariant record = transval.calculateVariant("PIK3R1:p.D464_R465del");
        //YDRL: TAT GAT AGA TTA, G of the D is at 68_293_799
        checkSingleHotspot(record, "TGATAGA", "T", "chr5", 68_293_798);
    }

    @Test
    public void delRangeTest()
    {
        TransvalVariant record = transval.calculateVariant("PIK3R1:p.D464_Y467del");
        //E   Y   D   R   L   Y   E
        //GAA TAT GAT AGA TTA TAT GAA

        assertEquals("ENST00000521381", record.transcriptId()); // canonical
        checkSingleHotspot(record, "AATATGATAGATT", "A", "chr5", 68_293_794);
    }

    @Test
    public void delInFirstExon()
    {
        TransvalVariant del = transval.calculateVariant("VHL:p.A5del");
        //RAE: AGG GCG GAG, G of the A is at 10_141_860
        checkSingleHotspot(del, "GGGC", "G", "chr3", 10_141_858);
    }

    @Test
    public void delFirstAndLastCodonsOfExonCrossExonBounds()
    {
        /*
        ch3, +
          10,146,514
        G|GT CAC CTT TGG CTC TTC ....CCA G|TG
             H
             115
         */
        TransvalVariant del = transval.calculateVariant("VHL:p.H115del");
        checkHotspots(del,
                hotspot("TCAC", "T", "chr3", 10_146_515),
                hotspot("GTCA", "G", "chr3", 10_146_514));

        del = transval.calculateVariant("VHL:p.L116_L118del");
        // CAC CTT TGG CTC TTC
        checkSingleHotspot(del,"ACCTTTGGCT", "A", "chr3", 10_146_517);

        // TG TAT ACT CTG   VYTL
        del = transval.calculateVariant("VHL:p.T157del");
        checkSingleHotspot(del,"ATAC", "A", "chr3", 10_149_790);
    }

    @Test
    public void delReverseStrand()
    {
        TransvalVariant record = transval.calculateVariant("BRAF:p.V600_R603del");
        checkSingleHotspot(record, "ATCGAGATTTCAC", "A", "chr7", 140753325);
    }

    @Test
    public void delReverseStrandFirstExonStart()
    {
        /*
        This is on -ve chr7.
        Positive strand:
        140,924,689....
        |                 140,924,703
        |                 |
        GCT CAG CGC CGC CAT
        5   4   3   2   1
        S   L   A   A   M
         */
        TransvalVariant a3del = transval.calculateVariant("BRAF:p.A3del");
        checkSingleHotspot(a3del, "GCGC", "G", "chr7", 140_924_694);

        TransvalVariant l4del = transval.calculateVariant("BRAF:p.L4del");
        checkSingleHotspot(l4del, "TCAG", "T", "chr7", 140_924_691);

        TransvalVariant a3l4del = transval.calculateVariant("BRAF:p.A3_L4del");
        checkSingleHotspot(a3l4del, "TCAGCGC", "T", "chr7", 140_924_691);

        TransvalVariant a2l4del = transval.calculateVariant("BRAF:p.A2_L4del");
        checkSingleHotspot(a2l4del, "TCAGCGCCGC", "T", "chr7", 140_924_691);
    }

    @Test
    public void delReverseStrandSecondExonStart()
    {
        /*
        chr, -
        140,850,198
        |                  140,850,212
        |                  |
        TTT GAT ATT CCA CAC
        51  50  49  48  47
        K   I   N   W   V
         */
        TransvalVariant del = transval.calculateVariant("BRAF:p.V47del");
        checkSingleHotspot(del, "ACAC", "A", "chr7", 140_850_209);

        del = transval.calculateVariant("BRAF:p.W48del");
        checkSingleHotspot(del, "TCCA", "T", "chr7", 140_850_206);

        del = transval.calculateVariant("BRAF:p.W48_I50del");
        checkSingleHotspot(del, "TGATATTCCA", "T", "chr7", 140_850_200);
    }

    @Test
    public void delReverseStrandFirstCodonAcrossExonBounds()
    {
        /*
        chr7, -
        140,808,047
        |                   140,808,062
        |                   |
        AAT TGG TTT CTT CTC T|...
        208 207 206 205 204 203
        I   P   K   K   E   G
         */
        TransvalVariant del = transval.calculateVariant("BRAF:p.E204del");
        checkSingleHotspot(del, "TCTC", "T", "chr7", 140_808_058);

        del = transval.calculateVariant("BRAF:p.E204_P207del");
        checkSingleHotspot(del, "TTGGTTTCTTCTC", "T", "chr7", 140_808_049);
    }

    @Test
    public void delReverseStrandFirstAndLastCodonsAcrossExonBounds()
    {
        /*
        chr7, -
        140,800,469
        |                 140,800,481
        |                |
        GAC AAA CAG CAA A..
        291             287
        V   F   L   L   D
         */
        TransvalVariant del = transval.calculateVariant("BRAF:p.F290del");
        checkSingleHotspot(del, "CAAA", "C", "chr7", 140_800_471);

        del = transval.calculateVariant("BRAF:p.L289del");
        checkSingleHotspot(del, "ACAG", "A", "chr7", 140_800_474);

        del = transval.calculateVariant("BRAF:p.L288_F290del");
        checkSingleHotspot(del, "CAAACAGCAA", "C", "chr7", 140_800_471);
    }

    @Test
    public void delReverseStrandFirstExonEnd()
    {
        /*
        This is on -ve chr7.
        Positive strand:
        140,924,566....
        |                   140,924,580
        |                  /
        CTC CTC CGG AAT GGC
        E   E   P   I   A
        46              42
         */
        TransvalVariant variant = transval.calculateVariant("BRAF:p.A42del");
        checkSingleHotspot(variant, "TGGC", "T", "chr7", 140_924_577);

        variant = transval.calculateVariant("BRAF:p.P44del");
        checkSingleHotspot(variant, "CCGG", "C", "chr7", 140_924_571);

        variant = transval.calculateVariant("BRAF:p.A42_P44del");
        checkSingleHotspot(variant, "TCCGGAATGG", "T", "chr7", 140_924_570);
        //
        //        TransvalVariant variant =  transval.calculateVariant("BRAF:p.E45del");
        //        checkSingleHotspot(variant,"CTC", "", "chr7", 140_924_566);
        //
        //        variant =  transval.calculateVariant("BRAF:p.E46del");
        //        checkSingleHotspot(variant,"CTC", "", "chr7", 140_924_566);
        //
        //        variant =  transval.calculateVariant("BRAF:p.E45_E46del");
        //        checkSingleHotspot(variant,"CTCCTC", "", "chr7", 140_924_566);
    }

    @Test
    public void deleteRangeInRegionOfRepeatingAminoAcids()
    {
        TransvalVariant variant = transval.calculateVariant("ARID1A:p.A345_A349del");
        checkHotspots(variant,
                hotspot("GGGCTGCGGCGGCGGC", "G", "chr1", 26697416),
                hotspot("GGCTGCGGCGGCGGCA", "G", "chr1", 26697417),
                hotspot("AGCTGCGGCGGCGGCC", "A", "chr1", 26697432),
                hotspot("TGCGGCGGCGGCCGCC", "T", "chr1", 26697435)
        );
    }

    @Test
    public void duplicationInRegionOfRepeatingAminoAcids()
    {
        TransvalVariant variant = transval.calculateVariant("ARID1A", "P20_P21dup");
        checkHotspots(variant,
                hotspot("A", "ACCCGCC", "chr1", 26_696_447),
                hotspot("C", "CCCGCCG", "chr1", 26_696_448),
                hotspot("G", "GCCGCCC", "chr1", 26_696_460)
        );
    }
}