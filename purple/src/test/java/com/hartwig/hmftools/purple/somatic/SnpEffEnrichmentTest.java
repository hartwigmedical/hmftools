package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SnpEffEnrichmentTest
{
    private final List<HmfTranscriptRegion> transcripts = HmfGenePanelSupplier.allGeneList37();
    private final Set<String> driverGenes = Sets.newHashSet();
    private VCFCodec codec;
    private SnpEffEnrichment victim;
    private List<VariantContext> capture;

    @Before
    public void setup()
    {
        codec = createTestCodec();
        capture = Lists.newArrayList();
        driverGenes.add("RP11-307N16.6");
        driverGenes.add("SPATA13");
        driverGenes.add("SPATA19");
        driverGenes.add("CUX1");
        driverGenes.add("CYLD");
        victim = new SnpEffEnrichment(driverGenes, transcripts, capture::add);
    }

    @NotNull
    private static VCFCodec createTestCodec()
    {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("SAMPLE"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @Test
    public void useFirstGeneIfNonInCanonical()
    {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;"
                        + "TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|"
                        + "nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|"
                        + "SPATA13|ENSG00000182957|transcript|ENST00000382108|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|"
                        + "3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final VariantImpact summary = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(0));
        assertEquals("RP11-307N16.6", summary.gene());
    }

    @Test
    public void favourCanonicalGeneWhenPossible()
    {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;"
                        + "TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|"
                        + "nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|"
                        + "LOW|SPATA13|ENSG00000182957|transcript|ENST00000424834|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|"
                        + "3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final VariantImpact summary = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(0));
        assertEquals("SPATA13", summary.gene());
    }

    @Test
    public void canonicalFieldsUseTranscriptAnnotation()
    {
        final String line =
                "11\t133715264\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;SOMATIC;TQSS=1;"
                        + "TQSS_NT=1;set=snvs;ANN=T|sequence_feature|MODERATE|SPATA19|ENSG00000166118|modified-residue:Phosphoserine|"
                        + "ENST00000299140|protein_coding|1/7|c.78G>A||||||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|"
                        + "ENSG00000166118|transcript|ENST00000299140|protein_coding|1/7|c.78G>A|p.Ser26Ser|133/861|78/504|26/167||,T|"
                        + "splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000532889|protein_coding|"
                        + "1/7|c.78G>A|p.Ser26Ser|170/653|78/504|26/167||\tGT:AD:DP\t0/1:57,49:108";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final VariantImpact summary = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(0));

        assertEquals("SPATA19", summary.gene());
        assertEquals(CodingEffect.SPLICE, summary.CanonicalCodingEffect);
        assertEquals("c.78G>A", summary.CanonicalHgvsCodingImpact);
        assertEquals("p.Ser26Ser", summary.CanonicalHgvsProteinImpact);
    }

    @Test
    public void testInframeIndelIgnoredForUTR()
    {
        final String template =
                "7\t101898665\t.\tCAA\tC\t767\tPASS\tANN=C|3_prime_UTR_variant|MODIFIER|CUX1|ENSG00000257923|transcript|ENST00000360264|"
                        + "protein_coding|24/24|c.*6344_*6345delAA|||||6344|,C|3_prime_UTR_variant|MODIFIER|CUX1|ENSG00000257923|"
                        + "transcript|ENST00000292535|protein_coding|24/24|c.*6344_*6345delAA|||||6344|,C|intron_variant|MODIFIER|CUX1|"
                        + "ENSG00000257923|transcript|ENST00000437600|protein_coding|14/22|c.1250-17971_1250-17970delAA||||||,C|"
                        + "intron_variant|MODIFIER|CUX1|ENSG00000257923|transcript|ENST00000292538|protein_coding|14/22|"
                        + "c.1256-17971_1256-17970delAA||||||,C|intron_variant|MODIFIER|CUX1|ENSG00000257923|transcript|ENST00000393824|"
                        + "protein_coding|13/21|c.1139-17971_1139-17970delAA||||||,C|intron_variant|MODIFIER|CUX1|ENSG00000257923|"
                        + "transcript|ENST00000547394|protein_coding|13/21|c.1208-17971_1208-17970delAA||||||,C|intron_variant|MODIFIER|"
                        + "CUX1|ENSG00000257923|transcript|ENST00000425244|protein_coding|13/21|c.1118-17971_1118-17970delAA||||||,C|"
                        + "intron_variant|MODIFIER|CUX1|ENSG00000257923|transcript|ENST00000560541|processed_transcript|14/22|"
                        + "n.1844-17971_1844-17970delAA||||||,C|intron_variant|MODIFIER|CUX1|ENSG00000257923|transcript|ENST00000558836|"
                        + "processed_transcript|13/21|n.1362-17971_1362-17970delAA||||||;LPS=13958;MAPPABILITY=1;MH=AA;%s"
                        + "RC=CCCAAAAAAAAAAAGAC;RC_MH=AA;RC_NM=2;RC_REPC=11;RC_REPS=A;REP_C=13;REP_S=A;SEC=CUX1,ENST00000360264,"
                        + "UTR_variant,MISSENSE,c.*6344_*6345delAA,;SEW=CUX1,ENST00000360264,UTR_variant,MISSENSE,1;TIER=LOW_CONFIDENCE;"
                        + "TNC=CCA\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:5,0:0:16:235,0:7,0:0,0,0,0,5,16:0:0,0,0:0,"
                        + "0,0,0,98,229:18\t0/1:14,70:0.515:136:1122,2161:32,69:29,9,32,0,14,136:0:0,1,2:628,141,593,0,259,2569:151";

        final String line1 = String.format(template, "");
        final String line2 = String.format(template, "PII=1;");

        final VariantContext variantContext1 = codec.decode(line1);
        final VariantContext variantContext2 = codec.decode(line2);
        victim.accept(variantContext1);
        victim.accept(variantContext2);
        final VariantImpact summary1 = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(0));
        final VariantImpact summary2 = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(1));
        assertEquals(CodingEffect.NONE, summary1.CanonicalCodingEffect);
        assertEquals(CodingEffect.NONE, summary2.CanonicalCodingEffect);
    }

    @Test
    public void testInframeIndel()
    {
        final String template =
                "16\t50813650\t.\tCCTCCTGTGAACTCACT\tC\t143\tPASS\tANN=C|frameshift_variant|HIGH|CYLD|ENSG00000083799|transcript|"
                        + "ENST00000311559|protein_coding|10/20|c.1214_1229delCTCCTGTGAACTCACT|p.Pro405fs|1605/5371|1214/2871|405/956||,"
                        + "C|frameshift_variant|HIGH|CYLD|ENSG00000083799|transcript|ENST00000569418|protein_coding|8/18|"
                        + "c.1205_1220delCTCCTGTGAACTCACT|p.Pro402fs|1483/3513|1205/2862|402/953||,C|frameshift_variant|HIGH|CYLD|"
                        + "ENSG00000083799|transcript|ENST00000540145|protein_coding|9/19|c.1214_1229delCTCCTGTGAACTCACT|p.Pro405fs|"
                        + "1629/8713|1214/2871|405/956||,C|frameshift_variant|HIGH|CYLD|ENSG00000083799|transcript|ENST00000564326|"
                        + "protein_coding|7/17|c.1205_1220delCTCCTGTGAACTCACT|p.Pro402fs|1399/3259|1205/2862|402/953||,C|"
                        + "frameshift_variant|HIGH|CYLD|ENSG00000083799|transcript|ENST00000566206|protein_coding|8/18|"
                        + "c.1205_1220delCTCCTGTGAACTCACT|p.Pro402fs|1455/3104|1205/2733|402/910||,C|frameshift_variant|HIGH|CYLD|"
                        + "ENSG00000083799|transcript|ENST00000398568|protein_coding|8/18|c.1205_1220delCTCCTGTGAACTCACT|p.Pro402fs|"
                        + "1505/3536|1205/2862|402/953||,C|frameshift_variant|HIGH|CYLD|ENSG00000083799|transcript|ENST00000427738|"
                        + "protein_coding|8/18|c.1214_1229delCTCCTGTGAACTCACT|p.Pro405fs|1419/8503|1214/2871|405/956||,C|sequence_feature|"
                        + "MODERATE|CYLD|ENSG00000083799|modified-residue:phosphoserine|ENST00000311559|protein_coding|10/20|"
                        + "c.1214_1229delCTCCTGTGAACTCACT||||||,C|sequence_feature|MODERATE|CYLD|ENSG00000083799|modified-residue:"
                        + "phosphoserine|ENST00000398568|protein_coding|8/18|c.1205_1220delCTCCTGTGAACTCACT||||||,C|intron_variant|"
                        + "MODIFIER|CYLD|ENSG00000083799|transcript|ENST00000568704|protein_coding|5/13|"
                        + "c.1129+1799_1129+1814delCTCCTGTGAACTCACT||||||,C|non_coding_transcript_exon_variant|MODIFIER|CYLD|"
                        + "ENSG00000083799|transcript|ENST00000569891|retained_intron|9/12|n.1600_1615delCTCCTGTGAACTCACT||||||,C|"
                        + "non_coding_transcript_exon_variant|MODIFIER|CYLD|ENSG00000083799|transcript|ENST00000563629|retained_intron|"
                        + "6/10|n.941_956delCTCCTGTGAACTCACT||||||;LOF=(CYLD|ENSG00000083799|17|0.41);LPS=6910;MAPPABILITY=1;MH;%s"
                        + "PURPLE_AF=0.0871;PURPLE_CN=3.97;PURPLE_GERMLINE=AMPLIFICATION;PURPLE_MACN=0.00;PURPLE_VCN=0.346;RC=ACGGCGACCAC;"
                        + "RC_NM=2;REPORTED;REP_C=3;REP_S=CCT;SEC=CYLD,ENST00000427738,frameshift_variant,MISSENSE,"
                        + "c.1214_1229delCTCCTGTGAACTCACT,p.Pro405fs;SEW=CYLD,ENST00000311559,frameshift_variant,MISSENSE,1;SUBCL=0.140;"
                        + "TIER=PANEL;TNC=TCC\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:36,1:0.026:38:1389,0:39,0:1,0,"
                        + "0,0,36,38:0:0,0,0:28,0,0,0,880,908:40\t0/1:119,6:0.047:128:4472,193:125,5:6,0,0,0,119,128:0:0,0,0:143,0,0,0,"
                        + "2582,2738:133";

        final String line1 = String.format(template, "");
        final String line2 = String.format(template, "PII=1;");

        final VariantContext variantContext1 = codec.decode(line1);
        final VariantContext variantContext2 = codec.decode(line2);
        victim.accept(variantContext1);
        victim.accept(variantContext2);
        final VariantImpact summary1 = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(0));
        final VariantImpact summary2 = SnpEffUtils.fromSnpEffEnrichedVariant(capture.get(1));
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, summary1.CanonicalCodingEffect);
        assertEquals(CodingEffect.MISSENSE, summary2.CanonicalCodingEffect);
    }
}
