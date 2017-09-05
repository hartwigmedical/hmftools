package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.VariantAnnotationTest.createVariantAnnotationBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.ImmutableHmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.ImmutableVariantAnnotation;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConsequenceDeterminerTest {

    private static final String CHROMOSOME = "X";
    private static final long POSITION = 150L;
    private static final long WRONG_POSITION = 1500L;
    private static final String GENE = "GENE";

    private static final String TRANSCRIPT = "TRANS";
    private static final int TRANSCRIPT_VERSION = 1;

    private static final String REF = "R";
    private static final String ALT = "A";
    private static final String COSMIC_ID = "123";
    private static final int ALLELE_READ_COUNT = 1;
    private static final int TOTAL_READ_COUNT = 2;

    private static final String HGVS_CODING = "c.RtoA";
    private static final String HGVS_PROTEIN = "p.RtoA";

    private static final String CHROMOSOME_BAND = "p1";
    private static final String ENTREZ_ID = "11";
    private static final String GENE_ID = "ENSG0000";
    private static final long GENE_START = 1;
    private static final long GENE_END = 42;

    @Test
    public void worksAsExpected() {
        final Slicer slicer = SlicerFactory.fromSingleGenomeRegion(
                ImmutableBEDGenomeRegion.of(CHROMOSOME, POSITION - 10, POSITION + 10));
        final Map<String, HmfGenomeRegion> transcriptMap = Maps.newHashMap();
        transcriptMap.put(TRANSCRIPT, hmfRegion());

        final ConsequenceDeterminer determiner = new ConsequenceDeterminer(slicer, transcriptMap);

        final VariantConsequence rightConsequence = VariantConsequence.MISSENSE_VARIANT;
        final VariantConsequence wrongConsequence = VariantConsequence.OTHER;

        final ImmutableVariantAnnotation.Builder annotationBuilder = createVariantAnnotationBuilder().featureID(TRANSCRIPT).
                featureType(ConsequenceDeterminer.FEATURE_TYPE_TRANSCRIPT).gene(GENE).hgvsCoding(HGVS_CODING).
                hgvsProtein(HGVS_PROTEIN);
        final VariantAnnotation rightAnnotation = annotationBuilder.consequences(
                Lists.newArrayList(rightConsequence)).build();
        final VariantAnnotation wrongAnnotation = annotationBuilder.consequences(
                Lists.newArrayList(wrongConsequence)).build();

        final SomaticVariant.Builder variantBuilder = new SomaticVariant.Builder().
                chromosome(CHROMOSOME).ref(REF).alt(ALT).cosmicID(COSMIC_ID).
                totalReadCount(TOTAL_READ_COUNT).alleleReadCount(ALLELE_READ_COUNT);

        final SomaticVariant rightVariant = variantBuilder.position(POSITION).
                annotations(Lists.newArrayList(rightAnnotation)).build();
        final SomaticVariant wrongConsequenceVariant = variantBuilder.position(POSITION).
                annotations(Lists.newArrayList(wrongAnnotation)).build();
        final SomaticVariant wrongPositionVariant = variantBuilder.position(WRONG_POSITION).
                annotations(Lists.newArrayList(rightAnnotation)).build();

        final List<VariantReport> findings = determiner.run(
                Lists.newArrayList(rightVariant, wrongConsequenceVariant, wrongPositionVariant)).findings();
        assertEquals(1, findings.size());

        final VariantReport report = findings.get(0);
        assertEquals(GENE, report.gene());
        assertEquals(CHROMOSOME + ":" + POSITION, report.chromosomePosition());
        assertEquals(REF, report.ref());
        assertEquals(ALT, report.alt());
        assertEquals(TRANSCRIPT + "." + TRANSCRIPT_VERSION, report.transcript());
        assertEquals(HGVS_CODING, report.hgvsCoding());
        assertEquals(HGVS_PROTEIN, report.hgvsProtein());
        assertEquals(rightConsequence.readableSequenceOntologyTerm(), report.consequence());
        assertEquals(COSMIC_ID, report.cosmicID());
        assertEquals(TOTAL_READ_COUNT, report.totalReadCount());
        assertEquals(ALLELE_READ_COUNT, report.alleleReadCount());
    }

    @NotNull
    private static HmfGenomeRegion hmfRegion() {
        return ImmutableHmfGenomeRegion.builder()
                .chromosome(CHROMOSOME)
                .start(POSITION - 10)
                .end(POSITION + 10)
                .gene(GENE)
                .transcriptID(TRANSCRIPT)
                .transcriptVersion(TRANSCRIPT_VERSION)
                .chromosomeBand(CHROMOSOME_BAND)
                .entrezId(ENTREZ_ID)
                .geneID(GENE_ID)
                .geneStart(GENE_START)
                .geneEnd(GENE_END)
                .build();
    }
}
