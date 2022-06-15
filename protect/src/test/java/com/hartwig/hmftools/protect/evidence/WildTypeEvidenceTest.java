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
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
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
import com.hartwig.hmftools.common.wildtype.WildTypeFactoryTest;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WildTypeEvidenceTest {

    @Test
    public void canDetermineWildType() {

        ReportableVariant variantGermline = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA1")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList(variantGermline);

        ReportableVariant variantSomatic = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA2")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList(variantSomatic);

        GainLoss reportableAmp = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
        GainLoss reportableDel = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = create("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        ReportableHomozygousDisruption homozygousDisruption = create("NRAS");
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        ReportableGeneDisruption geneDisruption = WildTypeFactoryTest.createDisruption("MYC");
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        //Test wild-type with somatic variant
        ActionableGene wildTypeSomaticVariant = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("BRCA2")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceSomaticVariant =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeSomaticVariant), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeSoamticVariant = wildTypeEvidenceSomaticVariant.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);

        assertEquals(evidencesWildTypeSoamticVariant.size(), 0);

        //Test wild-type with germline variant
        ActionableGene wildTypeGermlineVariant = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("BRCA1")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceGermlineVariant =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeGermlineVariant), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeGermlineVariant = wildTypeEvidenceGermlineVariant.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeGermlineVariant.size(), 0);

        //Test wild-type with CNV
        ActionableGene wildTypeCNV = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("APC")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceCNV =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeCNV), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeCNV = wildTypeEvidenceCNV.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeCNV.size(), 0);

        //Test wild-type with fusion  5 prime
        ActionableGene wildTypeFusion5 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("BAG4")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceFusion5 =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeFusion5), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeFusion5 = wildTypeEvidenceFusion5.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeFusion5.size(), 0);

        //Test wild-type with fusion  3 prime
        ActionableGene wildTypeFusion3 = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("BAG4")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceFusion3 =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeFusion3), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeFusion3 = wildTypeEvidenceFusion3.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeFusion3.size(), 0);

        //Test wild-type with homozygous disruption
        ActionableGene wildTypeHomozygousDisruption = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("NRAS")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceHomozygousDisruption =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeHomozygousDisruption), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeHomozygousDisruption = wildTypeEvidenceHomozygousDisruption.evidence(
                reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeHomozygousDisruption.size(), 0);

        //Test wild-type with gene disruption
        ActionableGene wildTypeGeneDisruption = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("MYC")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidenceGeneDisruption =
                new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildTypeGeneDisruption), listDriverGenes);

        List<ProtectEvidence> evidencesWildTypeGeneDisruption = wildTypeEvidenceGeneDisruption.evidence(
                reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildTypeGeneDisruption.size(), 0);

        //Test calling wild type
        ActionableGene wildType = ImmutableActionableGene.builder()
                .from(ServeTestFactory.createTestActionableGene())
                .gene("EGFR")
                .event(GeneLevelEvent.WILD_TYPE)
                .source(Knowledgebase.CKB)
                .build();

        WildTypeEvidence wildTypeEvidence = new WildTypeEvidence(EvidenceTestFactory.create(), Lists.newArrayList(wildType), listDriverGenes);

        List<ProtectEvidence> evidencesWildType = wildTypeEvidence.evidence(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions);
        assertEquals(evidencesWildType.size(), 1);
    }

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