package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Disruption;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Fusion;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.ImmutableFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableGeneFusionFactoryTest {

    private static final double EPSILON = 1.0e-10;

    public Fusion fusionTestData() {
        return ImmutableFusion.builder()
                .reportable(true)
                .knownType("3P-Prom")
                .primarySource("")
                .clusterId("")
                .clusterCount("")
                .resolvedType("")
                .svIdUp("9104868")
                .chrUp("7")
                .posUp("106696332")
                .orientUp("1")
                .typeUp("DEL")
                .ploidyUp(1.65)
                .geneUp("PRKAR2B")
                .chrBandUp("q22.3")
                .transcriptUp("ENST00000265717")
                .strandUp("1")
                .regionTypeUp("Intronic")
                .codingTypeUp("Coding")
                .exonUp("1")
                .phaseUp("1")
                .exonMaxUp("11")
                .disruptiveUp("true")
                .exactBaseUp("-1")
                .codingBasesUp("307")
                .totalCodingUp("1254")
                .codingStartUp("106685353")
                .codingEndUp("106800027")
                .transStartUp("106685094")
                .transEndUp("106802256")
                .distancePrevUp("10673")
                .canonicalUp("true")
                .biotypeUp("protein_coding")
                .svIdDown("15444932")
                .chrDown("7")
                .posDown("116411792")
                .orientDown("-1")
                .typeDown("DEL")
                .ploidyDown(1.65)
                .geneDown("MET")
                .chrBandDown("q31.2")
                .transcriptDown("ENST00000318493")
                .strandDown("1")
                .regionTypeDown("Intronic")
                .codingTypeDown("Coding")
                .exonDown("14")
                .phaseDown("1")
                .exonMaxDown("21")
                .disruptiveDown("true")
                .exactBaseDown("-1")
                .codingBasesDown("1283")
                .totalCodingDown("4224")
                .codingStartDown("116339139")
                .codingEndDown("116436178")
                .transStartDown("116312459")
                .transEndDown("116436396")
                .distancePrevDown("84")
                .canonicalDown("true")
                .biotypeDown("protein_coding")
                .proteinsKept("Protein kinase domain")
                .proteinsLost("Sema domain")
                .build();
    }

    @Test
    public void canConvertGeneFusions() {
      //  GeneFusion fusion1 =
       //         createFusion("TMPRSS2", "ENST00000398585", 4, 5, "PNPLA7", "ENST00000406427", 2, 3, KnownFusionsModel.CIVIC, 0.4);
       // GeneFusion fusion2 = createFusion("CLCN6", "ENST00000346436", 1, 2, "BRAF", "ENST00000288602", 8, 9, KnownFusionsModel.ONCOKB, 1D);

        List<ReportableGeneFusion> reportableFusions =
                ReportableGeneFusionFactory.fusionConvertToReportable(Lists.newArrayList(fusionTestData()));

        assertEquals(1, reportableFusions.size());

        ReportableGeneFusion reportableFusion1 = reportableFusions.get(0);
        assertEquals("PRKAR2B", reportableFusion1.geneStart());
        assertEquals("MET", reportableFusion1.geneEnd());
        assertEquals(1.65, reportableFusion1.ploidy(), EPSILON);
    }

    @NotNull
    private static GeneFusion createFusion(@NotNull String startGene, @NotNull String startTranscript, int startExonUpstream,
            int startExonDownstream, @NotNull String endGene, @NotNull String endTranscript, int endExonUpstream, int endExonDownstream,
            @NotNull String source, double ploidy) {
        Transcript upstreamTranscript = createFusionLeg(true, startGene, startTranscript, startExonUpstream, startExonDownstream, ploidy);
        Transcript downstreamTranscript = createFusionLeg(false, endGene, endTranscript, endExonUpstream, endExonDownstream, ploidy);

        return new GeneFusion(upstreamTranscript, downstreamTranscript, source, true);
    }

    @NotNull
    private static Transcript createFusionLeg(boolean isUpstream, @NotNull String gene, @NotNull String transcript, int exonUpstream,
            int exonDownstream, double ploidy) {
        EnrichedStructuralVariant variant = createEnrichedStructuralVariantBuilder().type(StructuralVariantType.BND)
                .primaryKey(1)
                .start(createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome("any").build())
                .end(createEnrichedStructuralVariantLegBuilder().orientation((byte) 1).chromosome("any").build())
                .ploidy(ploidy)
                .build();
        GeneAnnotation upstreamGene =
                new GeneAnnotation(variant, isUpstream, gene, Strings.EMPTY, 1, Lists.newArrayList(), Lists.newArrayList(), Strings.EMPTY);
        return new Transcript(upstreamGene, 0, transcript, exonUpstream, -1, exonDownstream, -1, 0, 0, 10, true, 0, 0, null, null);
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariantLeg.Builder createEnrichedStructuralVariantLegBuilder() {
        return ImmutableEnrichedStructuralVariantLeg.builder().homology(Strings.EMPTY).position(1);
    }

    @NotNull
    private static ImmutableEnrichedStructuralVariant.Builder createEnrichedStructuralVariantBuilder() {
        return ImmutableEnrichedStructuralVariant.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .qualityScore(0)
                .recovered(false);
    }
}