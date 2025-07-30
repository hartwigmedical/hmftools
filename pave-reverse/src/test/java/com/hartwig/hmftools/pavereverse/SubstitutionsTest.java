package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public final class SubstitutionsTest extends ReversePaveTestBase
{
    @Test
    public void mtorSNV()
    {
        // This example is based on TransvarConverterTest in the serve codebase.
        // The idea is to check that the current code agrees with Transval on
        // the key fields used by serve.
        /*
        05:47:12 - [DEBUG] - Converting transvar output line to TransvarRecord: 'MTOR:p.L2230V	ENST00000361445 (protein_coding)	MTOR	-	chr1:g.11122101A>C/c.6688T>G/p.L2230V	inside_[cds_in_exon_48]	CSQN=Missense;reference_codon=TTA;candidate_codons=GTA,GTC,GTG,GTT;candidate_mnv_variants=chr1:g.11122099_11122101delTAAinsGAC,chr1:g.11122099_11122101delTAAinsCAC,chr1:g.11122099_11122101delTAAinsAAC;aliases=ENSP00000354558;source=Ensembl'
05:47:12 - [DEBUG] - Interpreting transvar variant: 'TransvarRecord{transcript=ENST00000361445, chromosome=1, gdnaPosition=11122101, variantSpanMultipleExons=false, annotation=TransvarSnvMnv{gdnaRef=A, gdnaAlt=C, referenceCodon=TTA, candidateCodons=[GTA, GTC, GTG, GTT]}}'
05:47:12 - [DEBUG] - Converted 'MTOR|null|p.L2230V' to 4 hotspot(s)
05:47:12 - [INFO ] - Printing hotspots for 'MTOR:p.L2230V' on transcript null
05:47:12 - [INFO ] -  Hotspot{ref=A, alt=C, chromosome=chr1, position=11122101}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=GAC, chromosome=chr1, position=11122099}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=CAC, chromosome=chr1, position=11122099}
05:47:12 - [INFO ] -  Hotspot{ref=TAA, alt=AAC, chromosome=chr1, position=11122099}
         */
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("MTOR:p.L2230V");

        assertEquals("ENST00000361445", variant.transcriptName());
        assertEquals("1", variant.Chromosome);
        checkChanges(variant,
                basesChange("A", "C", "chr1", 11_122_101),
                basesChange("TAA", "GAC", "chr1", 11_122_099),
                basesChange("TAA", "CAC", "chr1", 11_122_099),
                basesChange("TAA", "AAC", "chr1", 11_122_099)
        );
    }

    @Test
    public void brafSNV()
    {
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF:p.V600E");
        assertEquals("7", variant.Chromosome);
        checkChanges(variant,
                basesChange("A", "T", "chr7", 140753336),
                basesChange("CA", "TT", "chr7", 140753335)
        );

        variant = reversePave.calculateProteinVariant("BRAF:p.F583W");
        checkSingleChange(variant, "AA", "CC", "chr7", 140_753_386);

        variant = reversePave.calculateProteinVariant("BRAF:p.I582W");
        checkSingleChange(variant, "TAT", "CCA", "chr7", 140_753_389);

        variant = reversePave.calculateProteinVariant("BRAF:p.I582V");
        checkChanges(variant,
                basesChange("TAT", "GAC", "chr7", 140_753_389),
                basesChange("TAT", "CAC", "chr7", 140_753_389),
                basesChange("TAT", "AAC", "chr7", 140_753_389),
                basesChange("T", "C", "chr7", 140_753_391)
        );

        variant = reversePave.calculateProteinVariant("BRAF:p.I582G");
        checkChanges(variant,
                basesChange("TAT", "ACC", "chr7", 140_753_389),
                basesChange("TAT", "CCC", "chr7", 140_753_389),
                basesChange("TAT", "GCC", "chr7", 140_753_389),
                basesChange("AT", "CC", "chr7", 140_753_390)
        );
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
        // exon1/exon2...SYRG/GH... AGG TAC CGA G/GT CAC  End of Exon1 is G, 10142187.
        // R is {CG*, AGG, AGA}. Need one of these from G**, *GT. Option is C -> G at 10142187
        BaseSequenceVariants record = reversePave.calculateProteinVariant("VHL:p.G114R");
        assertEquals("ENST00000256474", record.transcriptName());
        assertEquals("3", record.Chromosome);
        checkSingleChange(record, "G", "C", "chr3", 10_142_187);
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
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("TET2:p.Y1294A");
        assertEquals("ENST00000380013", variant.transcriptName()); // TransvarConvertTest has ENST00000540549
        assertEquals("4", variant.Chromosome);
        //        assertEquals(10_142_187, record.Position); // serve example has 10_183_871, which is from v37, I think

        checkChanges(variant,
                basesChange("TA", "GC", "chr4", 105_259_695),
                basesChange("TAC", "GCA", "chr4", 105_259_695),
                basesChange("TAC", "GCG", "chr4", 105_259_695),
                basesChange("TAC", "GCT", "chr4", 105_259_695)
        );
    }

    @Test
    public void snvAcrossExonPositiveStrand()
    {
        // ...N E E/E R T...  ...AAT GAA GA/G AGA ACT...  exon6/exon7, G at start of exon7 is at 105,259,619
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("TET2:p.E1268D");
        checkChanges(variant,
                basesChange("G", "T", "chr4", 105_259_619),
                basesChange("G", "C", "chr4", 105_259_619)
        );

        // ... A at the end of exon6 is at 105,243,778
        variant = reversePave.calculateProteinVariant("TET2:p.E1268V");
        checkChanges(variant,
                basesChange("A", "T", "chr4", 105_243_778)
        );

        variant = reversePave.calculateProteinVariant("TET2:p.E1268Q");
        checkChanges(variant,
                basesChange("G", "C", "chr4", 105_243_777)
        );

        // ...C V E/E Q I I...   ...TGT GTA G/AG CAA ATT exon3/exon4, G at end of exon 3 is at
        variant = reversePave.calculateProteinVariant("TET2:p.E1137Q");
        checkChanges(variant,
                basesChange("G", "C", "chr4", 105_237_351)
        );

        // ... A at start of exon 4 is at 105,241,339
        variant = reversePave.calculateProteinVariant("TET2:p.E1137D");
        checkChanges(variant,
                basesChange("G", "T", "chr4", 105_241_340),
                basesChange("G", "C", "chr4", 105_241_340)
        );
    }

    @Test
    public void snvAcrossExonNegativeStrand()
    {
        // exon15/exon14 == ...F I N/N N S... ==  ->...AAA TAT AT/T ATT ACT... == <-...TTT ATA TA/A TAA TGA...
        // end of exon15 is 140,753,393
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF:p.N581L");
        // For L want {TTA, TTG, CT*}. We have [A**] and [*AT] as possible within-exon variants
        assertTrue(variant.changes().isEmpty());

        // For K we have two options from changing the last base of the (reverse strand) codon.
        // Transvar only reports A->C. Not sure why.
        variant = reversePave.calculateProteinVariant("BRAF:p.N581K");
        checkChanges(variant,
                basesChange("A", "C", "chr7", 140_753_392),
                basesChange("A", "T", "chr7", 140_753_392)
        );

        // For D we have a single option: changing the (reverse strand) T at the beginning of the codon to G.
        variant = reversePave.calculateProteinVariant("BRAF:p.N581D");
        checkSingleChange(variant, "T", "C", "chr7", 140_754_187);
    }
}