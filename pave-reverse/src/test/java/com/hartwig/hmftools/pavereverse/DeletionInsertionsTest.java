package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.junit.Test;

public final class DeletionInsertionsTest extends ReversePaveTestBase
{
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
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("EGFR:p.L747_A750delinsP");
        assertEquals("ENST00000275493", record.transcriptName());
        assertEquals("7", record.Chromosome);
        //        assertEquals(55_174_776, record.Position);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(4, hotspots.size());
        assertTrue(hotspots.contains(basesChange("TTAAGAGAAGCA", "CCT", "chr7", 55174776)));
        assertTrue(hotspots.contains(basesChange("TTAAGAGAAGCA", "CCG", "chr7", 55174776)));
        assertTrue(hotspots.contains(basesChange("TTAAGAGAAG", "C", "chr7", 55174776)));
        assertTrue(hotspots.contains(basesChange("TTAAGAGAAGCA", "CCC", "chr7", 55174776)));
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
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("EGFR:p.A750_I759delinsG");
        assertEquals("ENST00000275493", record.transcriptName());
        assertEquals("7", record.Chromosome);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(4, hotspots.size());
        assertEquals(4, hotspots.size());
        assertTrue(hotspots.contains(basesChange("CAACATCTCCGAAAGCCAACAAGGAAATC", "GT", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(basesChange("CAACATCTCCGAAAGCCAACAAGGAAATC", "GG", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(basesChange("CAACATCTCCGAAAGCCAACAAGGAAATC", "GA", "chr7", 55_174_786)));
        assertTrue(hotspots.contains(basesChange("CAACATCTCCGAAAGCCAACAAGGAAAT", "G", "chr7", 55_174_786)));
    }

    @Test
    public void egfrServeExample2()
    {
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("EGFR:p.L747_K754delinsSPQ");
        assertEquals("ENST00000275493", record.transcriptName());
        assertEquals("7", record.Chromosome);
        // {TCA, AGC, AGT, TCC, TCG, TCT}*{CCA, CCC, CCG, CCT}*{CAA, CAG} we choose AGC CCA CAA, which bleeds into a section of As
        checkSingleChange(record, "TTAAGAGAAGCAACATCTCCGA", "AGCCCAC", "chr7", 55_174_776);
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
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("VHL:p.A5_W8delinsM");
        assertEquals("ENST00000256474", record.transcriptName());
        assertEquals("3", record.Chromosome);

        checkSingleChange(record, "GCGGAGAACTG", "AT", "chr3", 10141860);
    }

    @Test
    public void manyOptionsEGFR()
    {
        // EGFR exon 19 insertion EGFR-K745_E746insIPVAIK ...https://pmc.ncbi.nlm.nih.gov/articles/PMC10330422/
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("EGFR:p.K745_E746delinsIPVAIK");
        assertEquals("ENST00000275493", record.transcriptName());
        assertEquals("7", record.Chromosome);
        checkSingleChange(record, "AGG", "TACCAGTAGCAATAA", "chr7", 55174771);
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
        Because the inserted sequence has length > 2, we just choose the first alphabetically.

        This is something that Transvar doesn't handle:
        Converted 'ZYX|null|p.P67_D70delinsWKY' to 0 hotspot(s)
        Printing hotspots for 'ZYX:p.P67_D70delinsWKY' on transcript null
         */

        BaseSequenceVariants record =  reversePave.calculateProteinVariant("ZYX:p.P67_D70delinsWKY");
        assertEquals("ENST00000322764", record.transcriptName());
        assertEquals("7", record.Chromosome);
        checkSingleChange(record, "CCCCCGGAAG", "TGGAAAT", "chr7", 143_381_770);
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

        BaseSequenceVariants record =  reversePave.calculateProteinVariant("ZYX:p.E382_G385delinsDP");
        assertEquals("ENST00000322764", record.transcriptName());
        assertEquals("7", record.Chromosome);
        //        assertEquals(143_388_490, record.Position);
//        assertTrue(record.SpansMultipleExons);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(8, hotspots.size());
        // G|AA CTC TGC GGC -> G|AT CCC
        assertTrue(hotspots.contains(basesChange("ACTCTGCGG", "TCC", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "TCCT", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "TCCA", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "TCCG", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGG", "CCC", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "CCCT", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "CCCA", "chr7", 143_388_490)));
        assertTrue(hotspots.contains(basesChange("ACTCTGCGGC", "CCCG", "chr7", 143_388_490)));
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
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("VHL:p.R3_R4delinsQ");
        assertEquals("ENST00000256474", record.transcriptName());
        assertEquals("3", record.Chromosome);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(basesChange("GGAG", "A", "chr3", 10_141_855)));
        assertTrue(hotspots.contains(basesChange("GGAGG", "AA", "chr3", 10_141_855)));
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
        BaseSequenceVariants record =  reversePave.calculateProteinVariant("BRAF:p.A3_L4delinsE");
        assertEquals("ENST00000646891", record.transcriptName()); // canonical
        assertEquals("7", record.Chromosome);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(basesChange("AGCG", "T", "chr7", 140_924_693)));
        assertTrue(hotspots.contains(basesChange("CAGCG", "TT", "chr7", 140_924_692)));
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
        BaseSequenceVariants record = reversePave.calculateProteinVariant("PIK3R1:p.E451_Y452delinsD");
        assertEquals("ENST00000521381", record.transcriptName()); // canonical
        assertEquals("5", record.Chromosome);

        Set<BaseSequenceChange> hotspots = record.changes();
        assertEquals(2, hotspots.size());
        assertTrue(hotspots.contains(basesChange("AAT", "", "chr5", 68_293_761)));
        assertTrue(hotspots.contains(basesChange("ATAT", "C", "chr5", 68_293_762)));
    }

    @Test
    public void reverseStrandDelinsTest()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "N486_T491delinsK");
        checkChanges(variant,
                basesChange("GTAGGTGCTGTCACA", "", "chr7", 140778036),
                basesChange("TGTAGGTGCTGTCACA", "C", "chr7", 140778035)
        );
    }

    @Test
    public void delinsReverseStrandAtExonEndTest()
    {
        // -> ...TGC CAC ATC AC|C (...<- ACG GTG TAG TG|G ...A V D G|G)
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "G478_V480delinsE");
        // so -> GTG TAG TG|G goes to {AA|G, GA|G}. So delete TG TAG TG
        checkChanges(variant,
                basesChange("ACATCAC", "T", "chr7", 140778069),
                basesChange("CACATCAC", "TT", "chr7", 140778068)
        );
    }

    @Test
    public void delinsReverseStrandAtExonStartTest()
    {
        // -> ...C|CC AAT AGA GTC...CAA A|TC
        // <- ...G|GG TTA TCT CAG...GTT T|AG
        //       G    I   S   D...
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF", "D324_G327delinsQ");
        checkChanges(variant, basesChange("CCAATAGAGTC", "TG", "chr7", 140800362));
    }

    @Test
    public void length2Delins()
    {
        // Only the variants with the shortest edit distance should be returned.
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("EGFR", "L858_A859delinsRS");
        // CTG GCC -> {CG*, AGG, AGA}{TC*, AGC, AGT}. Variant with fewest changes is CGG TCC.
        Set<BaseSequenceChange> hotspots = variant.changes();
        assertEquals(1, hotspots.size());
        BaseSequenceChange hotspot = hotspots.iterator().next();
        assertEquals("TGG", hotspot.Ref);
        assertEquals("GGT", hotspot.Alt);
    }
}