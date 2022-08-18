package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SomaticVariantFactoryTest {

    private static final String SAMPLE = "sample";
    private static final String SOMATIC_VARIANT_FILE = Resources.getResource("variant/somatics.vcf").getPath();

    private static final double EPSILON = 1.0e-10;

    private SomaticVariantFactory victim;
    private VCFCodec codec;

    @Before
    public void setup() {
        victim = new SomaticVariantFactory();
        codec = createTestCodec();
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @Test
    public void canLoadSomaticSimpleVCFFromFile() throws IOException {
        final List<SomaticVariant> unfiltered = new SomaticVariantFactory().fromVCFFile("sample", SOMATIC_VARIANT_FILE);
        assertEquals(3, unfiltered.size());

        final List<SomaticVariant> filtered = SomaticVariantFactory.passOnlyInstance().fromVCFFile("sample", SOMATIC_VARIANT_FILE);
        assertEquals(2, filtered.size());
    }

    @Test
    public void canReadCorrectSomaticVariant() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:59,60:120";

        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("15", variant.chromosome());
        assertEquals(12345678, variant.position());
        assertEquals(VariantType.SNP, variant.type());

        assertEquals("C", variant.ref());
        assertEquals("A,G", variant.alt());

        assertEquals(0.5, variant.alleleFrequency(), EPSILON);
        assertEquals(120, variant.totalReadCount(), EPSILON);
    }

    @Test
    public void handleZeroReadCount() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:0,0:0";
        final Optional<SomaticVariant> missingAFVariant = victim.createVariant(SAMPLE, codec.decode(line));
        assertFalse(missingAFVariant.isPresent());
    }

    @Test
    public void incorrectSampleFieldYieldsMissingReadCounts() {
        final String missingAFLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:DP\t0/1:21";

        final Optional<SomaticVariant> missingAFVariant = victim.createVariant(SAMPLE, codec.decode(missingAFLine));
        assertFalse(missingAFVariant.isPresent());

        final String missingRefCovLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60:121";
        final Optional<SomaticVariant> missingRefCovVariant = victim.createVariant(SAMPLE, codec.decode(missingRefCovLine));
        assertFalse(missingRefCovVariant.isPresent());
    }

    @Test
    public void favourCanonicalGeneWhenPossible() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tSEC=SPATA13,ENST00000424834,synonymous_variant,SYNONYMOUS,c.3441C>T,p.Leu1147Leu;"
                        + "SEW=RP11-307N16.6,ENST00000382141,synonymous_variant,SYNONYMOUS,2;AC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;"
                        + "QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|"
                        + "ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|"
                        + "3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|transcript|ENST00000424834|"
                        + "protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("SPATA13", variant.gene());
    }

    /*
    @Test
    public void useFirstGeneIfNonInCanonical() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tSEW=RP11-307N16.6,ENST00000382141,synonymous_variant,SYNONYMOUS,2;AC=0;AF=0;AN=0;"
                        + "MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|"
                        + "synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|"
                        + "c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|"
                        + "transcript|ENST00000382108|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||"
                        + "\tGT:AD:DP\t0/1:36,38:75";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("RP11-307N16.6", variant.gene());
    }
     */

    @Test
    public void canonicalFieldsUseTranscriptAnnotation() {
        final String line =
                "11\t133715264\t.\tC\tT\t.\tPASS\tSEC=SPATA19,ENST00000299140,splice_region_variant&synonymous_variant,SYNONYMOUS,"
                        + "c.78G>A,p.Ser26Ser;SEW=SPATA19,ENST00000299140,splice_region_variant&synonymous_variant,SYNONYMOUS,1;AC=0;"
                        + "AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|"
                        + "sequence_feature|MODERATE|SPATA19|ENSG00000166118|modified-residue:Phosphoserine|ENST00000299140|protein_coding"
                        + "|1/7|c.78G>A||||||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|"
                        + "ENST00000299140|protein_coding|1/7|c.78G>A|p.Ser26Ser|133/861|78/504|26/167||,T|splice_region_variant&"
                        + "synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000532889|protein_coding|1/7|c.78G>A|"
                        + "p.Ser26Ser|170/653|78/504|26/167||\tGT:AD:DP\t0/1:57,49:108";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("SPATA19", variant.gene());
        assertEquals(CodingEffect.SYNONYMOUS, variant.canonicalCodingEffect());
        assertEquals("c.78G>A", variant.canonicalHgvsCodingImpact());
        assertEquals("p.Ser26Ser", variant.canonicalHgvsProteinImpact());
    }

    @NotNull
    private static SomaticVariant assertedGet(@NotNull Optional<SomaticVariant> variant) {
        assert variant.isPresent();
        return variant.get();
    }
}