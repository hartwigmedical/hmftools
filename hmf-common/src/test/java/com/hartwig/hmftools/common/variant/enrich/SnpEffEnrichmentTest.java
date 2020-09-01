package com.hartwig.hmftools.common.variant.enrich;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummary;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffSummaryFactory;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SnpEffEnrichmentTest {

    private final List<CanonicalTranscript> transcripts = CanonicalTranscriptFactory.create37();
    private final DriverGenePanel genePanel = DriverGenePanelFactoryTest.testGenePanel();
    private VCFCodec codec;
    private SnpEffEnrichment victim;
    private List<VariantContext> capture;

    @Before
    public void setup() {
        codec = createTestCodec();
        capture = Lists.newArrayList();
        victim = new SnpEffEnrichment(genePanel.driverGenes(), transcripts, capture::add);
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("SAMPLE"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @Test
    public void useFirstGeneIfNonInCanonical() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|transcript|ENST00000382108|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final SnpEffSummary summary = SnpEffSummaryFactory.fromSage(capture.get(0));
        assertEquals("RP11-307N16.6", summary.gene());
    }

    @Test
    public void favourCanonicalGeneWhenPossible() {
        final String line =
                "13\t24871731\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=43;QSS_NT=43;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|synonymous_variant|LOW|RP11-307N16.6|ENSG00000273167|transcript|ENST00000382141|nonsense_mediated_decay|12/16|c.3075C>T|p.Leu1025Leu|3653/4157|3075/3318|1025/1105||,T|synonymous_variant|LOW|SPATA13|ENSG00000182957|transcript|ENST00000424834|protein_coding|11/13|c.3441C>T|p.Leu1147Leu|3763/8457|3441/3834|1147/1277||\tGT:AD:DP\t0/1:36,38:75";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final SnpEffSummary summary = SnpEffSummaryFactory.fromSage(capture.get(0));
        assertEquals("SPATA13", summary.gene());
    }

    @Test
    public void canonicalFieldsUseTranscriptAnnotation() {
        final String line =
                "11\t133715264\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;MAPPABILITY=1.000000;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs;ANN=T|sequence_feature|MODERATE|SPATA19|ENSG00000166118|modified-residue:Phosphoserine|ENST00000299140|protein_coding|1/7|c.78G>A||||||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000299140|protein_coding|1/7|c.78G>A|p.Ser26Ser|133/861|78/504|26/167||,T|splice_region_variant&synonymous_variant|LOW|SPATA19|ENSG00000166118|transcript|ENST00000532889|protein_coding|1/7|c.78G>A|p.Ser26Ser|170/653|78/504|26/167||\tGT:AD:DP\t0/1:57,49:108";
        final VariantContext variantContext = codec.decode(line);
        victim.accept(variantContext);
        final SnpEffSummary summary = SnpEffSummaryFactory.fromSage(capture.get(0));

        assertEquals("SPATA19", summary.gene());
        assertEquals(CodingEffect.SYNONYMOUS, summary.canonicalCodingEffect());
        assertEquals("c.78G>A", summary.canonicalHgvsCodingImpact());
        assertEquals("p.Ser26Ser", summary.canonicalHgvsProteinImpact());
    }

}
