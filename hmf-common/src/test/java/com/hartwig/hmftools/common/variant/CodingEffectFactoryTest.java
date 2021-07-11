package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.PROTEIN_PROTEIN_CONTACT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_ACCEPTOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STRUCTURAL_INTERACTION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;

import static junit.framework.TestCase.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class CodingEffectFactoryTest {

    private final CodingEffectFactory victim = new CodingEffectFactory(HmfGenePanelSupplier.allGeneList37());
    private final VariantContext dummyVariant =
            VariantContextFromString.decode("17\t7579312\t.\tC\tA\t2453\tPASS\tTIER=PANEL;TNC=CCG\tGT\t0/0\t0/1");

    @Test
    public void testDonorMinusOne() {
        String vcf = "17\t7579312\t.\tC\tA\t2453\tPASS\tANN=A|splice_region_variant&synonymous_variant|LOW|TP53|ENSG00000141510|transcript|"
                + "ENST00000269305|protein_coding|4/11|c.375G>T|p.Thr125Thr|565/2579|375/1182|125/393||,;TIER=PANEL;TNC=CCG\t"
                + "GT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:34,0:0:34:1384,0:36,0:0,0,0,0,34,34:0:0,0,0:0,0,0,0,"
                + "835,835:36\t0/1:1,97:0.96:101:41,3754:1,100:78,12,7,0,1,101:0:0,0,0:2300,153,62,0,32,2560:102";
        assertEffect(CodingEffect.SPLICE, vcf);
    }

    @Test
    public void testDonorMinusOneMnv() {
        String vcf = "17\t7579311\t.\tCC\tAA\t2453\tPASS\tANN=A|splice_region_variant&synonymous_variant|LOW|TP53|ENSG00000141510|"
                + "transcript|ENST00000269305|protein_coding|4/11|c.375G>T|p.Thr125Thr|565/2579|375/1182|125/393||,;TIER=PANEL;TNC=CCG\t"
                + "GT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:34,0:0:34:1384,0:36,0:0,0,0,0,34,34:0:0,0,0:0,0,0,0,"
                + "835,835:36\t0/1:1,97:0.96:101:41,3754:1,100:78,12,7,0,1,101:0:0,0,0:2300,153,62,0,32,2560:102";
        assertEffect(CodingEffect.SPLICE, vcf);
    }

    @Test
    public void testDonorPlusFive() {
        String vcf = "4\t74021343\t.\tA\tT\t824\tPASS\tANN=T|splice_region_variant&intron_variant|LOW|ANKRD17|ENSG00000132466|transcript|"
                + "ENST00000358602|protein_coding|5/33|c.1000+5T>A||||||,T|splice_region_variant&intron_variant|LOW|ANKRD17|"
                + "ENSG00000132466|transcript|ENST00000558247|protein_coding|5/33|c.652+5T>A||||||WARNING_TRANSCRIPT_NO_START_CODON,T|"
                + "splice_region_variant&intron_variant|LOW|ANKRD17|ENSG00000132466|transcript|ENST00000509867|protein_coding|5/33|"
                + "c.661+5T>A||||||,T|splice_region_variant&intron_variant|LOW|ANKRD17|ENSG00000132466|transcript|ENST00000330838|"
                + "protein_coding|5/32|c.1000+5T>A||||||,T|splice_region_variant&intron_variant|LOW|ANKRD17|ENSG00000132466|transcript|"
                + "ENST00000560372|nonsense_mediated_decay|3/22|c.451+5T>A||||||WARNING_TRANSCRIPT_NO_START_CODON,T|intron_variant|"
                + "MODIFIER|ANKRD17|ENSG00000132466|transcript|ENST00000561029|protein_coding|1/12|c.385-13834T>A||||||"
                + "WARNING_TRANSCRIPT_INCOMPLETE;MAPPABILITY=1;PURPLE_AF=0.3288;PURPLE_CN=3.09;PURPLE_GERMLINE=DIPLOID;PURPLE_MACN=1.04;"
                + "PURPLE_VCN=1.02;RC=GCTTG;RC_NM=1;SEC=ANKRD17,ENST00000358602,splice_donor_variant&splice_region_variant&intron_variant,"
                + "SPLICE,c.1000+5T>A,;SEW=ANKRD17,ENST00000358602,splice_donor_variant&splice_region_variant&intron_variant,SPLICE,1;"
                + "TIER=HIGH_CONFIDENCE;TNC=CAT\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:24,0:0:24:927,0:24,0:0,0,0,"
                + "0,24,24:0:0,0,0:0,0,0,0,668,668:24\t0/1:129,29:0.184:158:5051,1140:133,30:25,2,2,0,129,158:0:0,0,0:779,45,24,0,"
                + "3829,4677:163";
        assertEffect(CodingEffect.SPLICE, vcf);
    }

    @Test
    public void testAcceptorPlusThreeForwardStrandG() {
        String vcf = "5\t76640675\t.\tC\tG\t972\tPASS\tANN=G|splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|transcript|"
                + "ENST00000264917|protein_coding|6/21|c.798-3C>G||||||,G|splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|"
                + "transcript|ENST00000340978|protein_coding|6/20|c.798-3C>G||||||,G|splice_region_variant&intron_variant|LOW|PDE8B|"
                + "ENSG00000113231|transcript|ENST00000346042|protein_coding|6/18|c.798-3C>G||||||,G|splice_region_variant&intron_variant|"
                + "LOW|PDE8B|ENSG00000113231|transcript|ENST00000333194|protein_coding|6/20|c.798-3C>G||||||,G|splice_region_variant&"
                + "intron_variant|LOW|PDE8B|ENSG00000113231|transcript|ENST00000342343|protein_coding|5/20|c.738-3C>G||||||,G|"
                + "splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|transcript|ENST00000503963|protein_coding|4/6|c.84-3C>G"
                + "||||||WARNING_TRANSCRIPT_NO_STOP_CODON;BIALLELIC;MAPPABILITY=1;PURPLE_AF=1.0846;PURPLE_CN=1.06;PURPLE_GERMLINE=DIPLOID;"
                + "PURPLE_MACN=9.422e-03;PURPLE_VCN=1.15;RC=CTTTGAG;RC_NM=1;REP_C=3;REP_S=T;SEC=PDE8B,ENST00000264917,"
                + "splice_region_variant&intron_variant,NONE,c.798-3C>G,;SEW=PDE8B,ENST00000264917,splice_region_variant&intron_variant,"
                + "NONE,1;TIER=HIGH_CONFIDENCE;TNC=TCA\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:93,0:0:93:3851,0:96,"
                + "0:0,0,0,0,93,93:0:0,0,0:0,0,0,0,2757,2757:96\t0/1:8,35:0.814:43:331,1548:8,39:27,7,1,0,8,43:0:0,0,0:900,72,18,0,217,"
                + "1207:47";
        assertEffect(SPLICE, vcf);
    }

    @Test
    public void testAcceptorPlusThreeForwardStrandC() {
        String vcf = "5\t76640675\t.\tT\tC\t972\tPASS\tANN=G|splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|transcript|"
                + "ENST00000264917|protein_coding|6/21|c.798-3C>G||||||,G|splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|"
                + "transcript|ENST00000340978|protein_coding|6/20|c.798-3C>G||||||,G|splice_region_variant&intron_variant|LOW|PDE8B|"
                + "ENSG00000113231|transcript|ENST00000346042|protein_coding|6/18|c.798-3C>G||||||,G|splice_region_variant&intron_variant|"
                + "LOW|PDE8B|ENSG00000113231|transcript|ENST00000333194|protein_coding|6/20|c.798-3C>G||||||,G|splice_region_variant&"
                + "intron_variant|LOW|PDE8B|ENSG00000113231|transcript|ENST00000342343|protein_coding|5/20|c.738-3C>G||||||,G|"
                + "splice_region_variant&intron_variant|LOW|PDE8B|ENSG00000113231|transcript|ENST00000503963|protein_coding|4/6|c.84-3C>G"
                + "||||||WARNING_TRANSCRIPT_NO_STOP_CODON;BIALLELIC;MAPPABILITY=1;PURPLE_AF=1.0846;PURPLE_CN=1.06;PURPLE_GERMLINE=DIPLOID;"
                + "PURPLE_MACN=9.422e-03;PURPLE_VCN=1.15;RC=CTTTGAG;RC_NM=1;REP_C=3;REP_S=T;SEC=PDE8B,ENST00000264917,"
                + "splice_region_variant&intron_variant,NONE,c.798-3C>G,;SEW=PDE8B,ENST00000264917,splice_region_variant&intron_variant,"
                + "NONE,1;TIER=HIGH_CONFIDENCE;TNC=TCA\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:93,0:0:93:3851,0:96,"
                + "0:0,0,0,0,93,93:0:0,0,0:0,0,0,0,2757,2757:96\t0/1:8,35:0.814:43:331,1548:8,39:27,7,1,0,8,43:0:0,0,0:900,72,18,0,217,"
                + "1207:47";
        assertEffect(NONE, vcf);
    }

    @Test
    public void testAcceptorPlusThreeReverseStrandG() {
        String vcf = "13\t28610183\t.\tA\tG\t409.0\tPASS\tRC=CTGTA;RC_IDX=2;RC_LF=TGAGGTTTCC;RC_NM=1;RC_RF=GAAAAGAACG;REP_C=2;REP_S=TA;"
                + "TIER=PANEL;TNC=TAT;AC=0;AN=0;MAPPABILITY=1;ANN=G|splice_region_variant&intron_variant|LOW|FLT3|ENSG00000122025|"
                + "transcript|ENST00000380982|protein_coding|10/23|c.1310-3T>C||||||,G|splice_region_variant&intron_variant|LOW|FLT3|"
                + "ENSG00000122025|transcript|ENST00000241453|protein_coding|10/23|c.1310-3T>C||||||,G|splice_region_variant&"
                + "intron_variant|LOW|FLT3|ENSG00000122025|transcript|ENST00000380987|nonsense_mediated_decay|10/24|c.1310-3T>C||||||,G|"
                + "splice_region_variant&intron_variant|LOW|FLT3|ENSG00000122025|transcript|ENST00000537084|protein_coding|10/22|"
                + "c.1310-3T>C||||||\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t./.:17,18:0.514:35:579,666:17,18:14,3,0,"
                + "0,1,17,35:0:0,0,0:382,27,0,0,1,423,832:35\t./.:50,65:0.565:115:1824,2490:50,68:50,12,3,0,0,50,115:1:0,0,0:1545,"
                + "189,65,0,0,1385,3184:118";
        assertEffect(NONE, vcf);
    }

    @Test
    public void testAcceptorPlusThreeReverseStrandC() {
        String vcf = "13\t28610183\t.\tA\tC\t409.0\tPASS\tRC=CTGTA;RC_IDX=2;RC_LF=TGAGGTTTCC;RC_NM=1;RC_RF=GAAAAGAACG;REP_C=2;REP_S=TA;"
                + "TIER=PANEL;TNC=TAT;AC=0;AN=0;MAPPABILITY=1;ANN=G|splice_region_variant&intron_variant|LOW|FLT3|ENSG00000122025|"
                + "transcript|ENST00000380982|protein_coding|10/23|c.1310-3T>C||||||,G|splice_region_variant&intron_variant|LOW|FLT3|"
                + "ENSG00000122025|transcript|ENST00000241453|protein_coding|10/23|c.1310-3T>C||||||,G|splice_region_variant&"
                + "intron_variant|LOW|FLT3|ENSG00000122025|transcript|ENST00000380987|nonsense_mediated_decay|10/24|c.1310-3T>C||||||,G|"
                + "splice_region_variant&intron_variant|LOW|FLT3|ENSG00000122025|transcript|ENST00000537084|protein_coding|10/22|"
                + "c.1310-3T>C||||||\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t./.:17,18:0.514:35:579,666:17,18:14,3,0,0,"
                + "1,17,35:0:0,0,0:382,27,0,0,1,423,832:35\t./.:50,65:0.565:115:1824,2490:50,68:50,12,3,0,0,50,115:1:0,0,0:1545,189,65,0,"
                + "0,1385,3184:118";
        assertEffect(SPLICE, vcf);
    }

    @Test
    public void testSingleEffect() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, FRAMESHIFT_VARIANT);
        assertEffect(NONSENSE_OR_FRAMESHIFT, STOP_GAINED);

        assertEffect(MISSENSE, MISSENSE_VARIANT);
        assertEffect(MISSENSE, PROTEIN_PROTEIN_CONTACT);
        assertEffect(MISSENSE, STRUCTURAL_INTERACTION_VARIANT);
        assertEffect(MISSENSE, INFRAME_DELETION);
        assertEffect(MISSENSE, INFRAME_INSERTION);

        assertEffect(SPLICE, SPLICE_ACCEPTOR_VARIANT);
        assertEffect(SPLICE, SPLICE_DONOR_VARIANT);
        assertEffect(SYNONYMOUS, SYNONYMOUS_VARIANT);

        assertEffect(NONE, SPLICE_REGION_VARIANT);
        assertEffect(NONE, INTRON_VARIANT);
    }

    @Test
    public void testEffectPriority() {
        assertEffect(NONSENSE_OR_FRAMESHIFT, STOP_GAINED, MISSENSE_VARIANT, SPLICE_ACCEPTOR_VARIANT, INTRON_VARIANT);
        assertEffect(SPLICE, MISSENSE_VARIANT, SPLICE_ACCEPTOR_VARIANT, SYNONYMOUS_VARIANT, INTRON_VARIANT);
        assertEffect(MISSENSE, MISSENSE_VARIANT, SYNONYMOUS_VARIANT, INTRON_VARIANT);
        assertEffect(SPLICE, SPLICE_ACCEPTOR_VARIANT, INTRON_VARIANT);
        assertEffect(SYNONYMOUS, SYNONYMOUS_VARIANT, INTRON_VARIANT);
        assertEffect(NONE, INTRON_VARIANT);
    }

    private void assertEffect(@NotNull CodingEffect expected, @NotNull String line) {
        VariantContext variant = VariantContextFromString.decode(line);
        SnpEffAnnotation snpEffSummary = SnpEffAnnotationParser.fromContext(variant).get(0);
        CodingEffect codingEffect = victim.effect(variant, snpEffSummary.gene(), snpEffSummary.consequences());
        assertEquals(expected, codingEffect);
    }

    private void assertEffect(@NotNull final CodingEffect expected, @NotNull final VariantConsequence... consequences) {
        assertEquals(expected, victim.effect(dummyVariant, "dummy_gene", Lists.newArrayList(consequences)));
    }
}
