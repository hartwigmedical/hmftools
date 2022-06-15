package com.hartwig.hmftools.protect.evidence;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class WildTypeEvidenceTest {

//    @Test
//    @Ignore
//    public void canDetermineWildType() {
//
//        ReportableVariant variantGermline = ImmutableReportableVariant.builder()
//                .from(ReportableVariantTestFactory.create())
//                .gene("BRCA1")
//                .chromosome("1")
//                .position(56412)
//                .ref("A")
//                .alt("C")
//                .build();
//        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList(variantGermline);
//
//        ReportableVariant variantSomatic = ImmutableReportableVariant.builder()
//                .from(ReportableVariantTestFactory.create())
//                .gene("BRCA2")
//                .chromosome("1")
//                .position(56412)
//                .ref("A")
//                .alt("C")
//                .build();
//        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList(variantSomatic);
//
//        GainLoss reportableAmp = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
//        GainLoss reportableDel = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
//        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);
//
//        LinxFusion reportedFusionMatch = create("BAG4", "FGFR1");
//        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
//
//        ReportableHomozygousDisruption homozygousDisruption = create("NRAS");
//        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);
//
//        List<DriverGene> listDriverGenes =
//                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR"));
//
//        //Test wild-type with somatic variant
//        ActionableGene wildTypeSomaticVariant = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("BRCA2")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidenceSomaticVariant =
//                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeSomaticVariant), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildTypeSoamticVariant = wildTypeEvidenceSomaticVariant.evidence(reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//
//        assertEquals(evidencesWildTypeSoamticVariant.size(), 1);
//
//        //Test wild-type with germline variant
//        ActionableGene wildTypeGermlineVariant = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("BRCA1")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidenceGermlineVariant =
//                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeGermlineVariant), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildTypeGermlineVariant = wildTypeEvidenceGermlineVariant.evidence(reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//        assertEquals(evidencesWildTypeGermlineVariant.size(), 0);
//
//        //Test wild-type with CNV
//        ActionableGene wildTypeCNV = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("APC")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidenceCNV =
//                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeCNV), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildTypeCNV = wildTypeEvidenceCNV.evidence(reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//        assertEquals(evidencesWildTypeCNV.size(), 0);
//
//        //Test wild-type with fusion
//        ActionableGene wildTypeFusion = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("FGFR1")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidenceFusion =
//                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeFusion), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildTypeFusion = wildTypeEvidenceFusion.evidence(reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//        assertEquals(evidencesWildTypeFusion.size(), 0);
//
//        //Test wild-type with homozygous disruption
//        ActionableGene wildTypeHomozygousDisruption = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("NRAS")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidenceHomozygousDisruption =
//                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeHomozygousDisruption), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildTypeHomozygousDisruption = wildTypeEvidenceHomozygousDisruption.evidence(
//                reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//        assertEquals(evidencesWildTypeHomozygousDisruption.size(), 1);
//
//        //Test calling wild type
//        ActionableGene wildType = ImmutableActionableGene.builder()
//                .from(ServeTestFactory.createTestActionableGene())
//                .gene("EGFR")
//                .event(GeneLevelEvent.WILD_TYPE)
//                .source(Knowledgebase.CKB)
//                .build();
//
//        WildTypeEvidence wildTypeEvidence = new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildType), listDriverGenes);
//
//        List<ProtectEvidence> evidencesWildType = wildTypeEvidence.evidence(reportableGermlineVariant,
//                reportableSomaticVariant,
//                reportableSomaticGainsLosses,
//                reportableFusions,
//                homozygousDisruptions);
//        assertEquals(evidencesWildType.size(), 1);
//    }

    @NotNull
    private static List<DriverGene> createDriverMap(@NotNull List<String> genes) {
        List<DriverGene> driverGeneList = Lists.newArrayList();
        for (String gene : genes) {
            driverGeneList.add(createDriverGene(gene));
        }
        return driverGeneList;
    }

    public static DriverGene createDriverGene(final String name) {
        return ImmutableDriverGene.builder()
                .gene(name)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .likelihoodType(TSG)
                .reportGermlineDisruption(true)
                .build();
    }

    @NotNull
    private static ReportableHomozygousDisruption create(@NotNull String gene) {
        return ImmutableReportableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    private static LinxFusion create(@NotNull String geneStart, @NotNull String geneEnd) {
        return linxFusionBuilder(geneStart, geneEnd).build();
    }

    @NotNull
    private static ImmutableLinxFusion.Builder linxFusionBuilder(@NotNull String geneStart, @NotNull String geneEnd) {
        return ImmutableLinxFusion.builder()
                .from(LinxTestFactory.createMinimalTestFusion())
                .geneStart(geneStart)
                .geneEnd(geneEnd)
                .reported(true);
    }
}