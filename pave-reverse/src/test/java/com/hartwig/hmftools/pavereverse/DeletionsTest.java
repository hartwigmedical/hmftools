package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public final class DeletionsTest extends ReversePaveTestBase
{
    @Test
    public void delSingleAminoAcidTest()
    {
        BaseSequenceVariants record = reversePave.calculateProteinVariant("PIK3R1:p.K459del");
        // EK: GAA AAA, G is at 68_293_781
        assertEquals("ENST00000521381", record.transcriptName()); // canonical
        assertEquals("5", record.Chromosome);
        checkSingleChange(record, "GAAA", "G", "chr5", 68_293_781);

        record = reversePave.calculateProteinVariant("PIK3R1:p.D464del");
        //YDR: TAT GAT AGA, G of the D is at 68_293_799
        checkChanges(record, basesChange("TATG", "T", "chr5", 68_293_796));
    }

    @Test
    public void delTwoAminoAcidsTest()
    {
        BaseSequenceVariants record = reversePave.calculateProteinVariant("PIK3R1:p.D464_R465del");
        //YDRL: TAT GAT AGA TTA, G of the D is at 68_293_799
        checkSingleChange(record, "TGATAGA", "T", "chr5", 68_293_798);
    }

    @Test
    public void delRangeTest()
    {
        BaseSequenceVariants record = reversePave.calculateProteinVariant("PIK3R1:p.D464_Y467del");
        //E   Y   D   R   L   Y   E
        //GAA TAT GAT AGA TTA TAT GAA

        assertEquals("ENST00000521381", record.transcriptName()); // canonical
        checkSingleChange(record, "AATATGATAGATT", "A", "chr5", 68_293_794);
    }

    @Test
    public void delInFirstExon()
    {
        BaseSequenceVariants del = reversePave.calculateProteinVariant("VHL:p.A5del");
        //RAE: AGG GCG GAG, G of the A is at 10_141_860
        checkSingleChange(del, "GGGC", "G", "chr3", 10_141_858);
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
        //        TransvalVariant del = transval.calculateVariant("VHL:p.H115del");
        //        checkHotspots(del,
        //                hotspot("TCAC", "T", "chr3", 10_146_515),
        //                hotspot("GTCA", "G", "chr3", 10_146_514));

        //        del = transval.calculateVariant("VHL:p.L116_L118del");
        // CAC CTT TGG CTC TTC
        //        checkSingleHotspot(del, "ACCTTTGGCT", "A", "chr3", 10_146_517);

        // TG TAT ACT CTG   VYTL
        var del = reversePave.calculateProteinVariant("VHL:p.T157del");
        checkSingleChange(del, "ATAC", "A", "chr3", 10_149_790);
    }

    @Test
    public void delReverseStrand()
    {
        BaseSequenceVariants record = reversePave.calculateProteinVariant("BRAF:p.V600_R603del");
        checkSingleChange(record, "ATCGAGATTTCAC", "A", "chr7", 140753325);
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
        BaseSequenceVariants a3del = reversePave.calculateProteinVariant("BRAF:p.A3del");
        checkSingleChange(a3del, "GCGC", "G", "chr7", 140_924_694);

        BaseSequenceVariants l4del = reversePave.calculateProteinVariant("BRAF:p.L4del");
        checkSingleChange(l4del, "TCAG", "T", "chr7", 140_924_691);

        BaseSequenceVariants a3l4del = reversePave.calculateProteinVariant("BRAF:p.A3_L4del");
        checkSingleChange(a3l4del, "TCAGCGC", "T", "chr7", 140_924_691);

        BaseSequenceVariants a2l4del = reversePave.calculateProteinVariant("BRAF:p.A2_L4del");
        checkSingleChange(a2l4del, "TCAGCGCCGC", "T", "chr7", 140_924_691);
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
        BaseSequenceVariants del = reversePave.calculateProteinVariant("BRAF:p.V47del");
        checkSingleChange(del, "ACAC", "A", "chr7", 140_850_209);

        del = reversePave.calculateProteinVariant("BRAF:p.W48del");
        checkSingleChange(del, "TCCA", "T", "chr7", 140_850_206);

        del = reversePave.calculateProteinVariant("BRAF:p.W48_I50del");
        checkSingleChange(del, "TGATATTCCA", "T", "chr7", 140_850_200);
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
        BaseSequenceVariants del = reversePave.calculateProteinVariant("BRAF:p.E204del");
        checkSingleChange(del, "TCTC", "T", "chr7", 140_808_058);

        del = reversePave.calculateProteinVariant("BRAF:p.E204_P207del");
        checkSingleChange(del, "TTGGTTTCTTCTC", "T", "chr7", 140_808_049);
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
        BaseSequenceVariants del = reversePave.calculateProteinVariant("BRAF:p.F290del");
        checkSingleChange(del, "CAAA", "C", "chr7", 140_800_471);

        del = reversePave.calculateProteinVariant("BRAF:p.L289del");
        checkSingleChange(del, "ACAG", "A", "chr7", 140_800_474);

        del = reversePave.calculateProteinVariant("BRAF:p.L288_F290del");
        checkSingleChange(del, "CAAACAGCAA", "C", "chr7", 140_800_471);
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
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("BRAF:p.A42del");
        checkSingleChange(variant, "TGGC", "T", "chr7", 140_924_577);

        variant = reversePave.calculateProteinVariant("BRAF:p.P44del");
        checkSingleChange(variant, "CCGG", "C", "chr7", 140_924_571);

        variant = reversePave.calculateProteinVariant("BRAF:p.A42_P44del");
        checkSingleChange(variant, "TCCGGAATGG", "T", "chr7", 140_924_570);
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
        BaseSequenceVariants variant = reversePave.calculateProteinVariant("ARID1A:p.A345_A349del");
        checkChanges(variant,
                basesChange("GGGCTGCGGCGGCGGC", "G", "chr1", 26697416),
                basesChange("GGCTGCGGCGGCGGCA", "G", "chr1", 26697417),
                basesChange("AGCTGCGGCGGCGGCC", "A", "chr1", 26697432),
                basesChange("TGCGGCGGCGGCCGCC", "T", "chr1", 26697435)
        );
    }
}