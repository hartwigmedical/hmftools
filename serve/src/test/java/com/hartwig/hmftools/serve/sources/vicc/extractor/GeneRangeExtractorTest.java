package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneRangeAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneRangeAnnotation;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class GeneRangeExtractorTest {

    @Test
    public void canExtractCodon() {
        assertEquals(Integer.valueOf("600"), GeneRangeExtractor.extractCodonNumber("BRAF (V600)"));
        assertEquals(Integer.valueOf("742"), GeneRangeExtractor.extractCodonNumber("W742"));
        assertEquals(Integer.valueOf("179"), GeneRangeExtractor.extractCodonNumber("Q179X"));
        assertEquals(Integer.valueOf("61"), GeneRangeExtractor.extractCodonNumber("KRAS Q61X"));
    }

    @Test
    public void canExtractExon() {
        assertEquals(Lists.newArrayList(19), GeneRangeExtractor.extractExonNumbers("EGFR exon 19 insertions"));
        assertEquals(Lists.newArrayList(20), GeneRangeExtractor.extractExonNumbers("ERBB2 proximal exon 20"));
        assertEquals(Lists.newArrayList(9, 11, 13, 14, 17), GeneRangeExtractor.extractExonNumbers("KIT mutation in exon 9,11,13,14 or 17"));
        assertEquals(Lists.newArrayList(16, 17, 18, 19), GeneRangeExtractor.extractExonNumbers("MET mutation in exon 16-19"));
        assertEquals(Lists.newArrayList(2, 3), GeneRangeExtractor.extractExonNumbers("Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(Lists.newArrayList(12), GeneRangeExtractor.extractExonNumbers("Exon 12 splice site insertion"));
    }

    //TODO
    @Test
    public void canExtractMutationFilter() {

    }

    //TODO
    @Test
    public void canExtractSpecificMutationTypeFilter(){

    }

    @Test
    public void canExtractRangesCodon() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("GNAQ", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "ENST00000286548",
                "GNAQ",
                "Q209",
                "description",
                "chromosome",
                "pos");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("GNAQ")
                        .chromosome("9")
                        .start(80409487)
                        .end(80409489)
                        .rangeInfo(209)
                        .mutationType(MutationTypeFilter.ANY)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    @Ignore
    public void canExtractRangesExon() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("GNAQ", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "ENST00000275493",
                "EGFR",
                "EXON 19 DELETION",
                "description",
                "chromosome",
                "pos");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("EGFR")
                        .chromosome("7")
                        .start(55242410)
                        .end(55242518)
                        .rangeInfo(19)
                        .exonId("ENSE00001756460")
                        .mutationType(MutationTypeFilter.MISSENSE_INFRAME_DELETION)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @Test
    public void canExtractRangesExonAndFusion() {
        GeneRangeExtractor geneRangeExtractor =
                new GeneRangeExtractor(HmfGenePanelSupplier.allGenesMap37(), createDriverGenes("GNAQ", "EGFR", "KIT"));
        Map<Feature, List<GeneRangeAnnotation>> geneRangesPerFeature = Maps.newHashMap();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "ENST00000288135",
                "KIT",
                "EXON 11 MUTATION",
                "description",
                "chromosome",
                "pos");

        geneRangesPerFeature.put(viccEntry.features().get(0),
                Lists.newArrayList(ImmutableGeneRangeAnnotation.builder()
                        .gene("KIT")
                        .chromosome("4")
                        .start(55593577)
                        .end(55593713)
                        .rangeInfo(11)
                        .exonId("ENSE00001074417")
                        .mutationType(MutationTypeFilter.MISSENSE_ANY)
                        .build()));

        assertEquals(geneRangesPerFeature, geneRangeExtractor.extractGeneRanges(viccEntry));
    }

    @NotNull
    private static List<DriverGene> createDriverGenes(@NotNull String gene1, @NotNull String gene2, @NotNull String gene3) {
        //TODO: determine real TSG and ONCO genes
        DriverGene driverGeneTsg = ImmutableDriverGene.builder()
                .gene(gene1)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(TSG)
                .build();

        DriverGene driverGeneOnco1 = ImmutableDriverGene.builder()
                .gene(gene2)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(ONCO)
                .build();

        DriverGene driverGeneOnco2 = ImmutableDriverGene.builder()
                .gene(gene3)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(ONCO)
                .build();
        return Lists.newArrayList(driverGeneTsg, driverGeneOnco1, driverGeneOnco2);
    }
}