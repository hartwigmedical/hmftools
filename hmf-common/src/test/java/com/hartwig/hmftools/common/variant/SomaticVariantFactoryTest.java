package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SomaticVariantFactoryTest {

    private static final String SAMPLE = "sample";
    private static final String BASE_PATH = Resources.getResource("variant").getPath();
    private static final String SOMATIC_VARIANT_FILE = "somatics.vcf";

    private static final double EPSILON = 1.0e-10;

    private SomaticVariantFactory victim;
    private VCFCodec codec;

    @Before
    public void setup() {
        victim = SomaticVariantFactory.unfilteredInstance();
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
        final String file = BASE_PATH + File.separator + SOMATIC_VARIANT_FILE;
        final List<SomaticVariant> unfiltered = SomaticVariantFactory.unfilteredInstance().fromVCFFile("sample", file);
        assertEquals(3, unfiltered.size());

        final List<SomaticVariant> filtered = SomaticVariantFactory.passOnlyInstance().fromVCFFile("sample", file);
        assertEquals(2, filtered.size());
    }

    @Test
    public void canReadCorrectSomaticVariant() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";

        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("15", variant.chromosome());
        assertEquals(12345678, variant.position());
        assertEquals(VariantType.SNP, variant.type());
        assertTrue(variant.isDBSNP());
        assertFalse(variant.isCOSMIC());

        assertEquals("C", variant.ref());
        assertEquals("A,G", variant.alt());

        assertEquals(0.5, variant.alleleFrequency(), EPSILON);
        assertEquals(120, variant.totalReadCount(), EPSILON);
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
    public void correctDBSNPAndCOSMIC() {
        final String both = "15\t12345678\trs1;COSM2\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasBoth = assertedGet(victim.createVariant(SAMPLE, codec.decode(both)));
        assertTrue(hasBoth.isDBSNP());
        assertTrue(hasBoth.isCOSMIC());
        assertTrue(hasBoth.cosmicIDs().contains("COSM2"));

        final String dbsnpOnly = "15\t12345678\trs1\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasDBSNPOnly = assertedGet(victim.createVariant(SAMPLE, codec.decode(dbsnpOnly)));
        assertTrue(hasDBSNPOnly.isDBSNP());
        assertFalse(hasDBSNPOnly.isCOSMIC());

        final String cosmicOnly = "15\t12345678\tCOSM2\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasCOSMICOnly = assertedGet(victim.createVariant(SAMPLE, codec.decode(cosmicOnly)));
        assertFalse(hasCOSMICOnly.isDBSNP());
        assertTrue(hasCOSMICOnly.isCOSMIC());
        assertTrue(hasCOSMICOnly.cosmicIDs().contains("COSM2"));

        final String none = "15\t12345678\tID\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasNone = assertedGet(victim.createVariant(SAMPLE, codec.decode(none)));
        assertFalse(hasNone.isDBSNP());
        assertFalse(hasNone.isCOSMIC());
    }

    @Test
    public void favourCanonicalGeneWhenPossible() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|transcript|ENST00000424834|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("SPATA13", variant.gene());
    }

    @Test
    public void useFirstGeneIfNonInCanonical() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|transcript|ENST00000382108|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("RP11-307N16.6", variant.gene());
    }

    @Test
    public void canonicalFieldsUseTranscriptAnnotation() {
        final String line =
                "11\t133715264\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|sequence_feature|MODERATE|SPATA19|ENSG00000166118|modified-residue:Phosphoserine|ENST00000299140|protein_coding|1/7|c.78G>A||||||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000299140|protein_coding|1/7|c.78G>A|p.Ser26Ser|133/861|78/504|26/167||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000532889|protein_coding|1/7|c.78G>A|p.Ser26Ser|170/653|78/504|26/167||\tGT:AD:DP\t0/1:57,49:108";
        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("SPATA19", variant.gene());
        assertEquals(CodingEffect.SYNONYMOUS, variant.canonicalCodingEffect());
        assertEquals("c.78G>A", variant.canonicalHgvsCodingImpact());
        assertEquals("p.Ser26Ser", variant.canonicalHgvsProteinImpact());
    }

    @Test
    public void revertToNormalCosmicIDIfCanonicalIsNoMatch() {
        final String sample = "sample";

        final String canonicalCosmicID =
                "14\t105246551\trs121434592;COSM33765\tC\tT\t.\tPASS\tCOSM2ENST=COSM33766|AKT1_ENST00000554581|c.49G>A|p.E17K|519\tGT:AD:DP\t0/1:120,71:204";

        final SomaticVariant canonicalCosmicIDVariant = assertedGet(SomaticVariantFactory.passOnlyInstance()
                .createVariant(sample, VariantContextFromString.decode(sample, canonicalCosmicID)));

        assertEquals("COSM33766", canonicalCosmicIDVariant.canonicalCosmicID());

        final String noCanonicalCosmicID =
                "14\t105246551\trs121434592;COSM33765\tC\tT\t.\tPASS\tCOSM2ENST=COSM33766|AKT1_ENST00000349310|c.49G>A|p.E17K|519\tGT:AD:DP\t0/1:120,71:204";

        final SomaticVariant noCanonicalCosmicIDVariant = assertedGet(SomaticVariantFactory.passOnlyInstance()
                .createVariant(sample, VariantContextFromString.decode(sample, noCanonicalCosmicID)));

        assertEquals("COSM33765", noCanonicalCosmicIDVariant.canonicalCosmicID());
    }

    @Test
    public void testFiltersAreReappliedAfterModifyingVariantContext() {
        victim = SomaticVariantFactory.unfilteredInstance();
        assertEquals(2, victim.process(SAMPLE, ponFilteredAndNearPonFiltered()).size());

        victim = SomaticVariantFactory.passOnlyInstance();
        assertEquals(0, victim.process(SAMPLE, ponFilteredAndNearPonFiltered()).size());
    }

    @NotNull
    private List<VariantContext> ponFilteredAndNearPonFiltered() {
        final VariantContext ponFiltered = codec.decode("15\t12345678\trs1;UCSC\tCAT\tA,G\t2\tGERMLINE_PON\tinfo;\tGT:AD:DP\t0/1:60,60:121");
        final VariantContext nearPonFiltered = codec.decode("15\t12345677\trs1;UCSC\tCAT\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60,60:121");
        return  Lists.newArrayList(nearPonFiltered, ponFiltered);
    }

    @NotNull
    private static SomaticVariant assertedGet(@NotNull Optional<SomaticVariant> variant) {
        assert variant.isPresent();
        return variant.get();
    }
}