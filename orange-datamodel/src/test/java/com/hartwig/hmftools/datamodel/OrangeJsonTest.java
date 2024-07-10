package com.hartwig.hmftools.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.time.LocalDate;
import java.util.Collection;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class OrangeJsonTest
{
    private static final String MINIMALLY_EMPTY_ORANGE_JSON =
            Thread.currentThread().getContextClassLoader().getResource("minimally.empty.orange.json").getPath();
    private static final String MINIMALLY_POPULATED_ORANGE_JSON =
            Thread.currentThread().getContextClassLoader().getResource("minimally.populated.orange.json").getPath();
    private static final String REAL_ORANGE_JSON = Thread.currentThread().getContextClassLoader().getResource("real.orange.json").getPath();

    private static final double EPSILON = 1.0E-2;

    @Test
    public void canReadMinimallyEmptyOrangeRecordJson() throws IOException
    {
        assertNotNull(OrangeJson.getInstance().read(MINIMALLY_EMPTY_ORANGE_JSON));
    }

    @Test
    public void canReadRealOrangeRecordJson() throws IOException
    {
        assertNotNull(OrangeJson.getInstance().read(REAL_ORANGE_JSON));
    }

    @Test
    public void canReadMinimallyPopulatedOrangeRecordJson() throws IOException
    {
        OrangeRecord record = OrangeJson.getInstance().read(MINIMALLY_POPULATED_ORANGE_JSON);

        assertEquals("TEST", record.sampleId());
        assertEquals(LocalDate.of(2022, 1, 20), record.samplingDate());
        assertEquals(OrangeRefGenomeVersion.V37, record.refGenomeVersion());

        assertPurple(record.purple());
        assertLinx(record.linx());
        assertPeach(record.peach());
        assertCuppa(record.cuppa());
        assertVirusInterpreter(record.virusInterpreter());
        assertLilac(record.lilac());
        assertChord(record.chord());
        assertPlots(record.plots());
    }

    private static void assertPurple(@NotNull PurpleRecord purple)
    {
        assertEquals(1, purple.fit().qc().status().size());
        assertTrue(purple.fit().qc().status().contains(PurpleQCStatus.PASS));
        assertEquals(0.12, purple.fit().purity(), EPSILON);
        assertEquals(3.1, purple.fit().ploidy(), EPSILON);

        assertEquals(0.1, purple.characteristics().microsatelliteIndelsPerMb(), EPSILON);
        assertEquals(PurpleMicrosatelliteStatus.MSS, purple.characteristics().microsatelliteStatus());
        assertEquals(13.71, purple.characteristics().tumorMutationalBurdenPerMb(), EPSILON);
        assertEquals(PurpleTumorMutationalStatus.HIGH, purple.characteristics().tumorMutationalBurdenStatus());
        assertEquals(185, purple.characteristics().tumorMutationalLoad());
        assertEquals(PurpleTumorMutationalStatus.HIGH, purple.characteristics().tumorMutationalLoadStatus());

        assertEquals(2, purple.somaticDrivers().size());
        PurpleDriver somaticDriver1 = findDriverByGene(purple.somaticDrivers(), "SF3B1");
        assertEquals("ENST00000335508", somaticDriver1.transcript());
        assertEquals(PurpleDriverType.MUTATION, somaticDriver1.type());
        assertEquals(0.2, somaticDriver1.driverLikelihood(), EPSILON);

        PurpleDriver somaticDriver2 = findDriverByGene(purple.somaticDrivers(), "SMAD4");
        assertEquals("ENST00000342988", somaticDriver2.transcript());
        assertEquals(PurpleDriverType.DEL, somaticDriver2.type());
        assertEquals(1.0, somaticDriver2.driverLikelihood(), EPSILON);

        assertEquals(1, purple.germlineDrivers().size());
        PurpleDriver germlineDriver1 = findDriverByGene(purple.germlineDrivers(), "BRCA1");
        assertEquals("ENST00000471181", germlineDriver1.transcript());
        assertEquals(PurpleDriverType.GERMLINE_MUTATION, germlineDriver1.type());
        assertEquals(0.8, germlineDriver1.driverLikelihood(), EPSILON);

        assertEquals(1, purple.allSomaticVariants().size());
        PurpleVariant somaticVariant = findVariantByGene(purple.allSomaticVariants(), "SF3B1");
        assertTrue(somaticVariant.reported());
        assertEquals(PurpleVariantType.SNP, somaticVariant.type());
        assertEquals("2", somaticVariant.chromosome());
        assertEquals(198266779, somaticVariant.position());
        assertEquals("G", somaticVariant.ref());
        assertEquals("A", somaticVariant.alt());
        assertEquals(2.03, somaticVariant.variantCopyNumber(), EPSILON);
        assertEquals(0.4, somaticVariant.minorAlleleCopyNumber(), EPSILON);
        assertEquals(3.02, somaticVariant.adjustedCopyNumber(), EPSILON);
        assertEquals(HotspotType.NON_HOTSPOT, somaticVariant.hotspot());
        assertEquals(20, somaticVariant.tumorDepth().totalReadCount());
        assertEquals(10, somaticVariant.tumorDepth().alleleReadCount());
        assertEquals(0.0, somaticVariant.subclonalLikelihood(), EPSILON);
        assertFalse(somaticVariant.biallelic());
        assertEquals(PurpleGenotypeStatus.UNKNOWN, somaticVariant.genotypeStatus());
        assertNull(somaticVariant.localPhaseSets());
        assertEquals("ENST00000335508", somaticVariant.canonicalImpact().transcript());
        assertEquals("c.2153C>T", somaticVariant.canonicalImpact().hgvsCodingImpact());
        assertEquals("p.Pro718Leu", somaticVariant.canonicalImpact().hgvsProteinImpact());
        assertEquals(2153, (int) somaticVariant.canonicalImpact().affectedCodon());
        assertEquals(12, (int) somaticVariant.canonicalImpact().affectedExon());
        assertFalse(somaticVariant.canonicalImpact().inSpliceRegion());
        assertEquals(1, somaticVariant.canonicalImpact().effects().size());
        assertTrue(somaticVariant.canonicalImpact().effects().contains(PurpleVariantEffect.MISSENSE));
        assertEquals(PurpleCodingEffect.MISSENSE, somaticVariant.canonicalImpact().codingEffect());
        assertTrue(somaticVariant.otherImpacts().isEmpty());

        assertEquals(1, purple.reportableSomaticVariants().size());
        assertEquals(somaticVariant, purple.reportableSomaticVariants().iterator().next());

        assertEquals(1, purple.allGermlineVariants().size());
        PurpleVariant germlineVariant = findVariantByGene(purple.allGermlineVariants(), "BRCA1");
        assertTrue(germlineVariant.reported());
        assertEquals(PurpleVariantType.SNP, germlineVariant.type());
        assertEquals("17", germlineVariant.chromosome());
        assertEquals(41209068, germlineVariant.position());
        assertEquals("C", germlineVariant.ref());
        assertEquals("T", germlineVariant.alt());
        assertEquals(1.0, germlineVariant.variantCopyNumber(), EPSILON);
        assertEquals(0.8, germlineVariant.minorAlleleCopyNumber(), EPSILON);
        assertEquals(2.0, germlineVariant.adjustedCopyNumber(), EPSILON);
        assertEquals(HotspotType.HOTSPOT, germlineVariant.hotspot());
        assertEquals(30, germlineVariant.tumorDepth().totalReadCount());
        assertEquals(20, germlineVariant.tumorDepth().alleleReadCount());
        assertEquals(0.2, germlineVariant.subclonalLikelihood(), EPSILON);
        assertFalse(germlineVariant.biallelic());
        assertEquals(PurpleGenotypeStatus.HET, germlineVariant.genotypeStatus());
        assertEquals(2, germlineVariant.localPhaseSets().size());
        assertTrue(germlineVariant.localPhaseSets().contains(1));
        assertTrue(germlineVariant.localPhaseSets().contains(2));
        assertEquals("ENST00000471181", germlineVariant.canonicalImpact().transcript());
        assertEquals("c.5340+1G>A", germlineVariant.canonicalImpact().hgvsCodingImpact());
        assertEquals("p.?", germlineVariant.canonicalImpact().hgvsProteinImpact());
        assertNull(germlineVariant.canonicalImpact().affectedCodon());
        assertNull(germlineVariant.canonicalImpact().affectedExon());
        assertTrue(germlineVariant.canonicalImpact().inSpliceRegion());
        assertEquals(2, germlineVariant.canonicalImpact().effects().size());
        assertTrue(germlineVariant.canonicalImpact().effects().contains(PurpleVariantEffect.SPLICE_DONOR));
        assertTrue(germlineVariant.canonicalImpact().effects().contains(PurpleVariantEffect.INTRONIC));
        assertEquals(PurpleCodingEffect.SPLICE, germlineVariant.canonicalImpact().codingEffect());
        assertTrue(germlineVariant.otherImpacts().isEmpty());

        assertEquals(1, purple.reportableGermlineVariants().size());
        assertEquals(germlineVariant, purple.reportableGermlineVariants().iterator().next());

        assertEquals(1, purple.allSomaticCopyNumbers().size());
        PurpleCopyNumber copyNumber = purple.allSomaticCopyNumbers().iterator().next();
        assertEquals("1", copyNumber.chromosome());
        assertEquals(10, copyNumber.start());
        assertEquals(20, copyNumber.end());
        assertEquals(4.1, copyNumber.averageTumorCopyNumber(), EPSILON);

        assertEquals(1, purple.allSomaticGeneCopyNumbers().size());
        PurpleGeneCopyNumber geneCopyNumber = purple.allSomaticGeneCopyNumbers().iterator().next();
        assertEquals("gene", geneCopyNumber.gene());
        assertEquals("12", geneCopyNumber.chromosome());
        assertEquals("p13", geneCopyNumber.chromosomeBand());
        assertEquals(1.2, geneCopyNumber.minCopyNumber(), EPSILON);
        assertEquals(0.4, geneCopyNumber.minMinorAlleleCopyNumber(), EPSILON);

        assertEquals(1, purple.allSomaticGainsLosses().size());
        PurpleGainLoss gainLoss = purple.allSomaticGainsLosses().iterator().next();
        assertEquals("5", gainLoss.chromosome());
        assertEquals("q2.2", gainLoss.chromosomeBand());
        assertEquals("SMAD4", gainLoss.gene());
        assertEquals("ENST00000591126", gainLoss.transcript());
        assertFalse(gainLoss.isCanonical());
        assertEquals(CopyNumberInterpretation.FULL_LOSS, gainLoss.interpretation());
        assertEquals(0.1, gainLoss.minCopies(), EPSILON);
        assertEquals(1.2, gainLoss.maxCopies(), EPSILON);

        assertEquals(1, purple.reportableSomaticGainsLosses().size());
        assertEquals(gainLoss, purple.reportableSomaticGainsLosses().iterator().next());
    }

    @NotNull
    private static PurpleDriver findDriverByGene(@NotNull Iterable<PurpleDriver> drivers, @NotNull String geneToFind)
    {
        for(PurpleDriver driver : drivers)
        {
            if(driver.gene().equals(geneToFind))
            {
                return driver;
            }
        }

        throw new IllegalStateException("Could not find driver for gene: " + geneToFind);
    }

    @NotNull
    private static PurpleVariant findVariantByGene(@NotNull Iterable<PurpleVariant> variants, @NotNull String geneToFind)
    {
        for(PurpleVariant variant : variants)
        {
            if(variant.gene().equals(geneToFind))
            {
                return variant;
            }
        }
        throw new IllegalStateException("Could not find variant for gene: " + geneToFind);
    }

    private static void assertLinx(@NotNull LinxRecord linx)
    {
        assertEquals(1, linx.allSomaticStructuralVariants().size());
        LinxSvAnnotation structuralVariant = linx.allSomaticStructuralVariants().iterator().next();
        assertEquals("id", structuralVariant.vcfId());
        assertEquals(1, structuralVariant.svId());
        assertEquals(2, structuralVariant.clusterId());
        assertEquals("", structuralVariant.clusterReason());
        assertFalse(structuralVariant.fragileSiteStart());
        assertFalse(structuralVariant.fragileSiteEnd());
        assertFalse(structuralVariant.isFoldback());
        assertEquals("NONE", structuralVariant.lineTypeStart());
        assertEquals("NONE", structuralVariant.lineTypeEnd());
        assertEquals(2, structuralVariant.junctionCopyNumberMin(), 0.01);
        assertEquals(3, structuralVariant.junctionCopyNumberMax(), 0.01);
        assertEquals("PTENR", structuralVariant.geneStart());
        assertEquals("PTEN", structuralVariant.geneEnd());
        assertEquals(0, structuralVariant.localTopologyIdStart());
        assertEquals(1, structuralVariant.localTopologyIdEnd());
        assertEquals("ISOLATED_S", structuralVariant.localTopologyStart());
        assertEquals("ISOLATED_BE", structuralVariant.localTopologyEnd());
        assertEquals(3, structuralVariant.localTICountStart());
        assertEquals(4, structuralVariant.localTICountEnd());

        assertEquals(1, linx.somaticHomozygousDisruptions().size());
        LinxHomozygousDisruption homozygousDisruption = linx.somaticHomozygousDisruptions().iterator().next();
        assertEquals("4", homozygousDisruption.chromosome());
        assertEquals("p1.12", homozygousDisruption.chromosomeBand());
        assertEquals("NF1", homozygousDisruption.gene());
        assertEquals("ENST00000358273", homozygousDisruption.transcript());
        assertTrue(homozygousDisruption.isCanonical());

        assertEquals(1, linx.allSomaticBreakends().size());
        LinxBreakend breakend = linx.allSomaticBreakends().iterator().next();
        assertFalse(breakend.reported());
        assertFalse(breakend.disruptive());
        assertEquals(1, breakend.svId());
        assertEquals("NF1", breakend.gene());
        assertEquals("1", breakend.chromosome());
        assertEquals("p12", breakend.chromosomeBand());
        assertEquals("trans", breakend.transcript());
        assertEquals(LinxBreakendType.DUP, breakend.type());
        assertEquals(1.1, breakend.junctionCopyNumber(), EPSILON);
        assertEquals(1.0, breakend.undisruptedCopyNumber(), EPSILON);
        assertEquals(-1, breakend.nextSpliceExonRank());
        assertEquals(1, breakend.exonUp());
        assertEquals(2, breakend.exonDown());
        assertEquals("Upstream", breakend.geneOrientation());
        assertEquals(-1, breakend.orientation());
        assertEquals(TranscriptRegionType.EXONIC, breakend.regionType());
        assertEquals(TranscriptCodingType.UTR_3P, breakend.codingType());

        assertEquals(1, linx.reportableSomaticBreakends().size());
        assertEquals(breakend, linx.reportableSomaticBreakends().iterator().next());

        assertEquals(1, linx.allSomaticFusions().size());
        LinxFusion fusion = linx.allSomaticFusions().iterator().next();
        assertTrue(fusion.reported());
        assertEquals(LinxFusionType.KNOWN_PAIR, fusion.reportedType());
        assertEquals(LinxUnreportableReason.NONE, fusion.unreportedReasons().iterator().next());
        assertEquals(1, fusion.unreportedReasons().size());
        assertEquals("TMPRSS2::ETV4", fusion.display());
        assertEquals("TMPRSS2", fusion.geneStart());
        assertEquals("ENST00000332149", fusion.geneTranscriptStart());
        assertEquals("Exon 1", fusion.geneContextStart());
        assertEquals(1, fusion.fusedExonUp());
        assertEquals("ETV4", fusion.geneEnd());
        assertEquals("ENST00000319349", fusion.geneTranscriptEnd());
        assertEquals("Exon 2", fusion.geneContextEnd());
        assertEquals(2, fusion.fusedExonDown());
        assertEquals(FusionLikelihoodType.HIGH, fusion.driverLikelihood());
        assertEquals(FusionPhasedType.INFRAME, fusion.phased());
        assertEquals(1.1, fusion.junctionCopyNumber(), EPSILON);

        assertEquals(1, linx.reportableSomaticFusions().size());
        assertEquals(fusion, linx.reportableSomaticFusions().iterator().next());
    }

    private static void assertPeach(@NotNull Collection<PeachGenotype> peach)
    {
        assertEquals(1, peach.size());
        PeachGenotype entry = peach.iterator().next();
        assertEquals("DPYD", entry.gene());
        assertEquals("*1", entry.allele());
        assertEquals(2, entry.alleleCount());
        assertEquals("*1_HOM", entry.haplotype());
        assertEquals("Normal Function", entry.function());
        assertEquals("5-Fluorouracil", entry.linkedDrugs());
        assertEquals("https://www.pharmgkb.org/guidelineAnnotation/PA166104939", entry.urlPrescriptionInfo());
        assertEquals("peach_prod_v1.3", entry.panelVersion());
        assertEquals("1.7", entry.repoVersion());
    }

    private static void assertCuppa(@NotNull CuppaData cuppa)
    {
        assertEquals(1, cuppa.predictions().size());
        CuppaPrediction prediction = cuppa.predictions().iterator().next();
        assertEquals("Melanoma", prediction.cancerType());
        assertEquals(0.996, prediction.likelihood(), EPSILON);
        assertEquals(0.979, prediction.snvPairwiseClassifier(), EPSILON);
        assertEquals(0.990, prediction.genomicPositionClassifier(), EPSILON);
        assertEquals(0.972, prediction.featureClassifier(), EPSILON);
        assertNull(prediction.altSjCohortClassifier());
        assertNull(prediction.expressionPairwiseClassifier());
    }

    private static void assertVirusInterpreter(@NotNull VirusInterpreterData virusInterpreter)
    {
        assertEquals(2, virusInterpreter.allViruses().size());
        VirusInterpreterEntry virus1 = findVirusByName(virusInterpreter.allViruses(), "Human papillomavirus 16");
        assertTrue(virus1.reported());
        assertEquals(VirusBreakendQCStatus.NO_ABNORMALITIES, virus1.qcStatus());
        assertEquals(VirusInterpretation.HPV, virus1.interpretation());
        assertEquals(1, virus1.integrations());
        assertEquals(VirusLikelihoodType.HIGH, virus1.driverLikelihood());
        assertEquals(0.9, virus1.percentageCovered(), EPSILON);

        VirusInterpreterEntry virus2 = findVirusByName(virusInterpreter.allViruses(), "Human betaherpesvirus 6B");
        assertFalse(virus2.reported());
        assertEquals(VirusBreakendQCStatus.NO_ABNORMALITIES, virus2.qcStatus());
        assertNull(virus2.interpretation());
        assertEquals(0, virus2.integrations());
        assertEquals(VirusLikelihoodType.LOW, virus2.driverLikelihood());
        assertEquals(0.4, virus2.percentageCovered(), EPSILON);

        assertEquals(1, virusInterpreter.reportableViruses().size());
        assertEquals(virus1, virusInterpreter.reportableViruses().iterator().next());
    }

    @NotNull
    private static VirusInterpreterEntry findVirusByName(@NotNull Iterable<VirusInterpreterEntry> entries, @NotNull String nameToFind)
    {
        for(VirusInterpreterEntry entry : entries)
        {
            if(entry.name().equals(nameToFind))
            {
                return entry;
            }
        }
        throw new IllegalStateException("Could not find virus with name: " + nameToFind);
    }

    private static void assertLilac(@NotNull LilacRecord lilac)
    {
        assertEquals("PASS", lilac.qc());

        assertEquals(1, lilac.alleles().size());
        LilacAllele allele = lilac.alleles().iterator().next();
        assertEquals("A*01:01", allele.allele());
        assertEquals(6.1, allele.tumorCopyNumber(), EPSILON);
        assertEquals(5.0, allele.somaticMissense(), EPSILON);
        assertEquals(4.0, allele.somaticNonsenseOrFrameshift(), EPSILON);
        assertEquals(3.0, allele.somaticSplice(), EPSILON);
        assertEquals(0.4, allele.somaticSynonymous(), EPSILON);
        assertEquals(1.0, allele.somaticInframeIndel(), EPSILON);
    }

    private static void assertChord(@NotNull ChordRecord chord)
    {
        assertEquals(0.01, chord.hrdValue(), EPSILON);
        assertEquals(ChordStatus.HR_PROFICIENT, chord.hrStatus());
    }

    private static void assertPlots(@NotNull OrangePlots plots)
    {
        assertEquals("plot/empty.circos.png", plots.purpleFinalCircosPlot());
    }
}