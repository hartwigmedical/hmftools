package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class HgvsCodingAddressesTest
{
    // all exons and introns of length 12
    String exon1 = "CTAGGACACGAG";
    String intron1 = "CGAGGGCCCAAA";
    String exon2 = "CGACACGAGTAA";
    String intron2 = "GGGCCCAAATTT";
    String exon3 = "ATGAAGGAACCT";
    String intron3 = "AAAGGGCCCTTT";
    String exon4 = "AAGATGGAACCT";
    String intron4 = "AACCGGTTACGT";
    String exon5 = "GTACACAAGCCT";
    String intron6 = "CACCGGTTACGT";
    String exon6 = "ATACACAAGCCT";

    private final String refBases = generateTestBases(20)
            + exon1 + intron1 + exon2 + intron2 + exon3 + intron3 + exon4 + intron4 + exon5 + intron6 + exon6
            + generateTestBases(20);
    MockRefGenome refGenome = new MockRefGenome();
    private final ImpactClassifier classifier;
    private final TranscriptData transcriptData;

    public HgvsCodingAddressesTest()
    {
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        classifier = new ImpactClassifier(refGenome);
        //            21           33           45           57           69           81           93           105          117          129          141
        //            |   exon1    |   intron1  |   exon2    |   intron2  |   exon3    |   intron4  |   exon4    |   intron5  |   exon5    |   intron6  |   exon6
        // random 20  CTAGGACACGAG CGAGGGCCCAAA CGACACGAGTAA GGGCCCAAATTT ATGAAGGAACCT AAAGGGCCCTTT AAGATGGAACCT AACCGGTTACGT GTACACAAGCCT CACCGGTTACGT ATACACAAGCCT random 20
        //            |                         |          |              |                         |                         |                         |
        //            -24                       -12        -1             1                         13                        *1                        *13
        transcriptData = createTransExons(GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 21, 45, 69, 93, 117, 141 }, 11, 69, 104, false, "");
    }

    @Test
    public void secondDownstreamIntronAndExon()
    {
        // boundary
//        checkDupAt(127, 2, "c.*12_*12+1dupTC");
        checkDelAt(127, 2, "c.*12_*12+1delTC");
        checkInsAt(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG

        // in intron6
        checkDupAt(128, 2, "c.*12+1_*12+2dupCA");
        checkDelAt(128, 2, "c.*12+1_*12+2delCA");
        checkInsAt(129, "GG",  "c.*12+1_*12+2insGG");
        checkSnvAt(129, "C", "G", "c.*12+1C>G");

        checkDupAt(133, 2, "c.*12+6_*12+7dupGT");
        checkDelAt(133, 2, "c.*12+6_*12+7delGT");
        checkInsAt(134, "AA",  "c.*12+6_*12+7insAA");
        checkSnvAt(134, "G", "A", "c.*12+6G>A");

//        checkDupAt(134, 2, "c.*13-6_*13-5dupTT");
//        checkDelAt(134, 2, "c.*13-6_*13-5delTT");
        checkInsAt(135, "AA",  "c.*13-6_*13-5insAA");
        checkSnvAt(135, "T", "A", "c.*13-6T>A");

        // exon6
        checkDupAt(140, 2, "c.*13_*14dupAT");
//        checkDelAt(140, 2, "c.*13_*14delAT");
        checkInsAt(141, "CC",  "c.*13_*14insCC");
        checkSnvAt(141, "A", "T", "c.*13A>T");

        checkDupAt(141, 2, "c.*14_*15dupTA");
        checkDelAt(141, 2, "c.*14_*15delTA");
        checkInsAt(142, "CC",  "c.*14_*15insCC");
        checkSnvAt(142, "T", "G", "c.*14T>G");

        checkDupAt(150, 2, "c.*23_*24dupCT");
        checkDelAt(150, 2, "c.*23_*24delCT");
        checkInsAt(151, "AA",  "c.*23_*24insAA");
        checkSnvAt(151, "C", "G", "c.*23C>G");

        checkSnvAt(152, "T", "G", "c.*24T>G");
    }

    @Test
    public void downstreamExon()
    {
        //        checkDupAt(116, 2, "c.*1_*2dupGT");
        //        checkDelAt(116, 2, "c.*1_*2delGT");
        checkInsAt(117, "CC", "c.*1_*2insCC");
        checkSnvAt(117, "G", "C", "c.*1G>C");

        checkDupAt(119, 2, "c.*4_*5dupCA");
        checkDelAt(119, 2, "c.*4_*5delCA");
        checkInsAt(120, "GG", "c.*4_*5insGG");
        checkSnvAt(120, "C", "A", "c.*4C>A");

        checkDupAt(126, 2, "c.*11_*12dupCT");
        checkDelAt(126, 2, "c.*11_*12delCT");
        checkInsAt(127, "GG", "c.*11_*12insGG");
        checkSnvAt(127, "C", "A", "c.*11C>A");

        checkDupAt(127, 1, "c.*12dupT");
        checkDelAt(127, 1, "c.*12delT");
        checkInsAt(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG ?
        checkSnvAt(128, "T", "A", "c.*12T>A");
    }

    @Test
    public void downnstreamIntronExonBoundary()
    {
        checkDupAt(115, 2, "c.*1-1_*1-0dupTG");
        //        checkDelAt(115, 2, "c.*1-1_*1delTG");
        //        checkInsAt(116, "CC", "c.*1-1_*1insCC");
    }

    @Test
    public void downstreamIntron()
    {
        checkDupAt(104, 2, "c.24+1_24+2dupAA");
        checkDelAt(104, 2, "c.24+1_24+2delAA");
        checkInsAt(105, "CC", "c.24+1_24+2insCC");
        checkSnvAt(105, "A", "C", "c.24+1A>C");

        checkDupAt(109, 2, "c.24+6_24+7dupGT"); //c.24+6_c.*1-6dupGT?
        checkDelAt(109, 2, "c.24+6_24+7delGT");
        checkInsAt(110, "CC", "c.24+6_24+7insCC");
        checkSnvAt(110, "G", "C", "c.24+6G>C");

        //        checkDupAt(109, 2, "c.*1-2_*1-1dupGT");
        checkDelAt(114, 2, "c.*1-2_*1-1delGT");
        checkInsAt(115, "CC", "c.*1-2_*1-1insCC");
        checkSnvAt(115, "G", "C", "c.*1-2G>C");
    }

    @Test
    public void secondExonDownstreamIntronBoundary()
    {
        //        checkDupAt(103, 2, "c.24_24+1dupTA");
        //        checkDelAt(103, 2, "c.24_24+1delTA");
        checkInsAt(104, "GG", "c.24-0_24+1insGG"); // c.24_24+1insGG?
    }

    @Test
    public void secondCodingExon()
    {
        checkDupAt(102, 2, "c.23_24dupCT");
        //        checkDelAt(102, 2, "c.23_24delCT");
        checkInsAt(103, "A", "c.23_24insA");
        checkSnvAt(103, "C", "G", "c.23C>G");

        checkDupAt(93, 3, "c.14_16dupAGA");
        checkDelAt(93, 3, "c.14_16delAGA");
        checkInsAt(94, "C", "c.14_15insC");
        checkSnvAt(94, "A", "C", "c.14A>C");

        //        checkDupAt(92, 3, "c.13_15dupAAG");
        checkDelAt(92, 3, "c.13_15delAAG");
        checkInsAt(93, "C", "c.13_14insC");
        checkSnvAt(93, "A", "C", "c.13A>C");
    }

    @Test
    public void firstIntronSecondExonBoundary()
    {
        checkDupAt(90, 3, "c.13-2_13-0dupTTA"); // c.13-2_13dupTTA ?
        checkDelAt(90, 3, "c.13-2_13delTTA");
    }

    @Test
    public void firstIntron()
    {
        checkDupAt(90, 2, "c.13-2_13-1dupTT");
        checkDelAt(90, 2, "c.13-2_13-1delTT");
        checkInsAt(91, "AC", "c.13-2_13-1insAC");
        checkSnvAt(92, "T", "A", "c.13-1T>A");

        //        checkDupAt(86, 3, "c.13-6_13-4dupCCC");
        //        checkDelAt(86, 3, "c.13-6_13-4delCCC");
        checkInsAt(87, "TT", "c.13-6_13-5insTT");
        checkSnvAt(87, "C", "T", "c.13-6C>T");

        checkDupAt(85, 3, "c.12+6_12+8dupGCC");
        checkDelAt(85, 3, "c.12+6_12+8delGCC");
        checkInsAt(86, "TT", "c.12+6_12+7insTT");
        checkSnvAt(86, "G", "T", "c.12+6G>T");

        checkDupAt(84, 3, "c.12+5_12+7dupGGC");
        checkDelAt(84, 3, "c.12+5_12+7delGGC");
        checkInsAt(85, "TT", "c.12+5_12+6insTT");
        checkSnvAt(85, "G", "T", "c.12+5G>T");

        checkDupAt(80, 2, "c.12+1_12+2dupAA");
        checkDelAt(80, 2, "c.12+1_12+2delAA");
        checkInsAt(81, "CC", "c.12+1_12+2insCC");
        checkSnvAt(81, "A", "C", "c.12+1A>C");
    }

    @Test
    public void firstExonFirstIntronBoundary()
    {
        //        checkDupAt(79, 2, "c.12_12+1dupTT");
        checkDelAt(79, 3, "c.12_12+2delTAA");
        checkInsAt(80, "CC", "c.12-0_12+1insCC"); // c.12_12+1insCC ??
    }

    @Test
    public void firstCodingExon()
    {
        checkDupAt(79, 1, "c.12dupT");
        checkDelAt(79, 1, "c.12delT");
        checkSnvAt(80, "T", "A", "c.12T>A");

        checkDupAt(78, 2, "c.11_12dupCT");
        checkDelAt(78, 2, "c.11_12delCT");
        checkInsAt(79, "AA", "c.11_12insAA");
        checkSnvAt(79, "C", "T", "c.11C>T");

        checkDupAt(76, 2, "c.9_10dupAC");
        checkDelAt(76, 2, "c.9_10delAC");
        checkInsAt(77, "TTT", "c.9_10insTTT");
        checkSnvAt(77, "A", "C", "c.9A>C");

        checkDupAt(71, 2, "c.4_5dupAA");
        checkDelAt(71, 2, "c.4_5delAA");
        checkInsAt(72, "TTT", "c.4_5insTTT");
        checkSnvAt(72, "A", "C", "c.4A>C");

        checkDupAt(69, 2, "c.2_3dupTG");
        checkDelAt(69, 2, "c.2_3delTG");
        checkInsAt(70, "C", "c.2_3insC");
        checkSnvAt(70, "T", "C", "c.2T>C");

        //        checkDupAt(68, 3, "c.1_3dupATG");
        checkDelAt(68, 3, "c.1_3delATG");
        checkInsAt(69, "CC", "c.1_2insCC");
        checkSnvAt(69, "A", "C", "c.1A>C");
        checkDelAt(68, 1, "c.1delA");
    }

    @Test
    public void firstUpstreamIntron()
    {
        //        checkDupAt(66, 2, "c.1-2_1-1dupTT");
        checkDelAt(66, 2, "c.1-2_1-1delTT");
        //        checkDupAt(64, 4, "c.1-4_1-1dupATTT");
        checkDelAt(64, 4, "c.1-4_1-1delATTT");
        //        checkDupAt(63, 3, "c.1-5_1-3dupAAT");
        checkDelAt(63, 3, "c.1-5_1-3delAAT");
        checkDelAt(63, 1, "c.1-5delA");
        checkDelAt(62, 1, "c.-1+7delA");
        checkDelAt(61, 1, "c.-1+6delC");
        checkDelAt(59, 1, "c.-1+4delC");
        checkDelAt(58, 1, "c.-1+3delG");
        checkDelAt(57, 1, "c.-1+2delG");
        checkDelAt(56, 1, "c.-1+1delG");
    }

    @Test
    public void firstUpstreamExon()
    {
        checkDelAt(55, 1, "c.-1delA");
        //        checkDupAt(62, 3, "c.-1+6_-1+8dupAAA"); // -1+6_1-4 ??
        //        checkDelAt(62, 3, "c.-1+6_-1+8delAAA");
        //        checkDupAt(59, 3, "c.-1+3_-1+5dupCCC");
        //        checkDelAt(59, 3, "c.-1+3_-1+5delCCC");
        //        checkDupAt(57, 3, "c.-1+1_-1+3dupGGC");
        //        checkDelAt(57, 3, "c.-1+1_-1+3delGGC");
        //        checkDupAt(56, 3, "c.-1+0_-1+2dupGGG");
        //        checkDelAt(56, 3, "c.-1_-1+2delGGG");
        //        checkDupAt(55, 2, "c.-2_-1dupAA");
        checkDelAt(55, 1, "c.-1delA");
        checkDelAt(55, 2, "c.-1_-1+1delAG");
        //        checkDupAt(54, 2, "c.-2_-1dupAA");
        checkDelAt(54, 2, "c.-2_-1delAA");
        //        checkDupAt(46, 2, "c.-11_-10dupAC");
        checkDelAt(46, 2, "c.-10_-9delAC");
        //        checkDupAt(45, 2, "c.-12_-11dupAC");
        checkDelAt(45, 2, "c.-11_-10delGA");
        checkDelAt(45, 1, "c.-11delG");
        //        checkDelAt(44, 1, "c.-12delC");
    }

    @Test
    public void secondUpstreamIntron()
    {
        checkDelAt(43, 1, "c.-12-1delA");
        checkDelAt(40, 4, "c.-12-4_-12-1delCAAA");
        checkDelAt(36, 2, "c.-13+5_-13+6delGG");
        checkDelAt(32, 2, "c.-13+1_-13+2delCG");
    }

    @Test
    public void secondUpstreamExon()
    {
        checkDelAt(31, 2, "c.-13_-13+1delGC");
        checkDelAt(27, 3, "c.-17_-15delACG");
        checkDelAt(23, 3, "c.-21_-19delGGA");
        checkDelAt(21, 1, "c.-23delT");
        //        checkDelAt(20, 1, "c.-24delC");
    }

    private void checkDupAt(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDupVariant(refBases, position, length));
    }

    private void checkSnvAt(int position, String ref, String alt, String expected)
    {
        checkBaseAt(position, ref);
        checkVariantImpact(expected, createSnvVariant(position, ref, alt));
    }

    private void checkBaseAt(int position, String ref)
    {
        assertEquals(refBases.substring(position - 1, position), ref);
    }

    private VariantData createDupVariant(String refBases, int position, int length)
    {
        int pos = position - 1;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + length + 1);
        VariantData var = new VariantData(CHR_1, position, ref, alt);
        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);
        return var;
    }

    private VariantData createSnvVariant(int position, String ref, String alt)
    {
        return new VariantData(CHR_1, position, ref, alt);

    }

    private void checkInsAt(int position, String bases, String expected)
    {
        checkVariantImpact(expected, createInsVariant(position, bases));
    }

    private void checkDelAt(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDelVariant(refBases, position, length));
    }

    private void checkVariantImpact(final String expected, final VariantData var)
    {
        VariantTransImpact impact = classifier.classifyVariant(var, transcriptData);
        assertEquals(expected, impact.codingContext().Hgvs);
    }

    private VariantData createInsVariant(int pos, String bases)
    {
        String ref = refBases.substring(pos, pos + 1);
        String alt = ref + bases;
        return new VariantData(CHR_1, pos, ref, alt);
    }

    private VariantData createDelVariant(String refBases, int position, int length)
    {
        int pos = position - 1;
        String ref = refBases.substring(pos, pos + length + 1);
        String alt = refBases.substring(pos, pos + 1);
        return new VariantData(CHR_1, position, ref, alt);
    }
}
