package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.TumorOnlyPurityAnalysis.snvTypeClassifier;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_38;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class TumorOnlyPurityAnalysisTest extends PurityTestBase
{
    private static final String TUMOR_ID = "tumor123";
    private static final int BASELINE_READINGS_PER_ARM = 100;
    private static final int BASELINE_DEPTH = 1000;
    private Map<Chromosome, SortedSet<AmberSite>> SitesByChromosome;
    private Map<HumanChromosome, SortedSet<PositionEvidence>> EvidenceByChromosome;
    private Map<PositionEvidence, Double> GnomadFrequencyByEvidence;
    private File OutputDir;
    private TumorOnlyPurityAnalysis Analysis;

    @Before
    public void setup() throws Exception
    {
        SitesByChromosome = new HashMap<>();
        EvidenceByChromosome = new HashMap<>();
        GnomadFrequencyByEvidence = new HashMap<>();
        OutputDir = Files.createTempDirectory("amber").toFile();
        OutputDir.deleteOnExit();
    }

    @Test
    public void mutationTypeClassifierTest()
    {
        assertEquals(CanonicalSnvType.C_A, snvTypeClassifier().apply(pe("C", "A", 0.5)));
        assertEquals(CanonicalSnvType.C_A, snvTypeClassifier().apply(pe("C", "A", 0.65)));
        assertEquals(CanonicalSnvType.T_G, snvTypeClassifier().apply(pe("C", "A", 0.66)));
    }

    @Test
    public void noPeaks()
    {
        createBaselineReadings();
        createAnalysis();
        assertEquals(0.04, Analysis.cutoff(), 0.001);
        assertTrue(Analysis.contaminationPeaks().isEmpty());
    }

    @Test
    public void onePeakAboveMinimumCutoff()
    {
        createBaselineReadings();
        final double vaf = 0.21;
        createCopyNumberEventPoints(vaf);
        createAnalysis();
        assertEquals(TumorOnlyPurityAnalysis.MIN_CUTOFF, Analysis.cutoff(), 0.001);
        assertEquals(0, Analysis.contaminationPeaks().size());
        assertEquals(1, Analysis.copyNumberPeaks().size());
        assertEquals(vaf, Analysis.copyNumberPeaks().get(0).vaf(), 0.01);
    }

    @Test
    public void oneCnEventBelowMinimumCutoff()
    {
        createBaselineReadings();
        final double vaf = 0.1;
        // Add a 20 point copy number peak on 1P.
        createCopyNumberEventPoints(vaf);
        createAnalysis();
        assertEquals(vaf / 3, Analysis.cutoff(), 0.01);
        assertEquals(0, Analysis.contaminationPeaks().size());
        assertEquals(2, Analysis.copyNumberPeaks().size());
        CandidatePeak peak0 = Analysis.copyNumberPeaks().get(0);
        assertEquals(1.0, peak0.homozygousProportion(), 0.01);
        assertEquals(vaf, peak0.vaf(), 0.01);
        assertEquals(20, peak0.homozygousEvidencePoints().size());
        assertEquals(0, peak0.heterozygousEvidencePoints().size());
        CandidatePeak peak1 = Analysis.copyNumberPeaks().get(1);
        assertEquals(vaf * 2.0, peak1.vaf(), 0.1);
        assertEquals(0.0, peak1.homozygousProportion(), 0.01);
        assertEquals(0, peak1.homozygousEvidencePoints().size());
        assertEquals(20, peak1.heterozygousEvidencePoints().size());
    }

    @Test
    public void contaminationTest()
    {
        createBaselineReadings();
        addInContaminationData(0.25);
        createAnalysis();
        assertEquals(1, Analysis.contaminationPeaks().size());
    }

    private void createCopyNumberEventPoints(final double vaf)
    {
        // Add a 20 point copy number peak on 1P.
        int position = startPosition(new ChrArm(HumanChromosome._1, Arm.P)) + 50; // offset so as not to interfere with baseline points
        boolean low = false;
        for(int i = 0; i < 20; i++)
        {
            addPeakPoint(HumanChromosome._1, position, BASELINE_DEPTH, vaf, 0.5, low);
            position += BASELINE_DEPTH;
            low = !low;
        }
    }

    private void createBaselineReadings()
    {
        for(ChrArm chrArm : chromosomeArms())
        {
            int position = startPosition(chrArm);
            for(int i = 0; i < BASELINE_READINGS_PER_ARM; i++)
            {
                addBaselinePoint(chrArm.chromosome(), position, BASELINE_DEPTH, 0.5);
                position += 1000;
            }
        }
    }

    private void addInContaminationData(double contamination)
    {
        int depth = (int) (BASELINE_DEPTH * contamination);
        for(HumanChromosome chromosome : EvidenceByChromosome.keySet())
        {
            for(PositionEvidence baselineEvidence : EvidenceByChromosome.get(chromosome))
            {
                double gnomadFrequency = GnomadFrequencyByEvidence.get(baselineEvidence);
                PositionEvidence contaminationEvidence =
                        createGermlinePositionEvidence(chromosome, baselineEvidence.position(), baselineEvidence.ref(), baselineEvidence.alt(), depth, gnomadFrequency);
                baselineEvidence.ReadDepth += contaminationEvidence.ReadDepth;
                baselineEvidence.RefSupport += contaminationEvidence.RefSupport;
                baselineEvidence.AltSupport += contaminationEvidence.AltSupport;

            }
        }
    }

    private void createAnalysis()
    {
        List<PositionEvidence> evidence = new ArrayList<>();
        for(HumanChromosome chromosome : EvidenceByChromosome.keySet())
        {
            evidence.addAll(EvidenceByChromosome.get(chromosome));
        }
        ListMultimap<Chromosome, AmberSite> sites = ArrayListMultimap.create();
        for(Chromosome chromosome : SitesByChromosome.keySet())
        {
            sites.putAll(chromosome, SitesByChromosome.get(chromosome));
        }
        PurityAnalysisConfig config = new PurityAnalysisConfig(TUMOR_ID, V38, OutputDir.getAbsolutePath(), 8);
        Analysis = new TumorOnlyPurityAnalysis(evidence, sites, config);
    }

    private int startPosition(ChrArm chrArm)
    {
        if(chrArm.arm() == Arm.P)
        {
            return 1_000_000;
        }
        int centromere = COORDS_38.centromeres().get(chrArm.chromosome());
        return (centromere + 3_000_000) / 1_000_000;
    }

    private List<ChrArm> chromosomeArms()
    {
        List<ChrArm> arms = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(!chromosome.hasShortArm())
            {
                arms.add(new ChrArm(chromosome, Arm.P));
            }
            arms.add(new ChrArm(chromosome, Arm.Q));
        }
        return arms;
    }

    private void addBaselinePoint(HumanChromosome chromosome, int position, int depth, double gnomadFrequencyOfSite)
    {
        Pair<String, String> refAndAlt = randomRefAndAltBases();
        String ref = refAndAlt.getLeft();
        String alt = refAndAlt.getRight();

        createAndAddAmberSite(chromosome, position, ref, alt, gnomadFrequencyOfSite);

        PositionEvidence evidence = createGermlinePositionEvidence(chromosome, position, ref, alt, depth, gnomadFrequencyOfSite);
        EvidenceByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(evidence);
        GnomadFrequencyByEvidence.put(evidence, gnomadFrequencyOfSite);
    }

    @NotNull
    private static PositionEvidence createGermlinePositionEvidence(final HumanChromosome chromosome, final int position,
            final String ref, final String alt,
            final int depth,
            final double gnomadFrequencyOfSite)
    {
        PositionEvidence evidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, ref, alt);
        int altDepth = altDepthForGermlinePoint(gnomadFrequencyOfSite, depth);
        evidence.ReadDepth = depth;
        evidence.RefSupport = depth - altDepth;
        evidence.AltSupport = altDepth;
        return evidence;
    }

    private void addPeakPoint(HumanChromosome chromosome, int position, int depth, double vaf, double gnomadFrequencyOfSite, boolean low)
    {
        Pair<String, String> refAndAlt = randomRefAndAltBases();
        String ref = refAndAlt.getLeft();
        String alt = refAndAlt.getRight();

        createAndAddAmberSite(chromosome, position, ref, alt, gnomadFrequencyOfSite);

        PositionEvidence evidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, ref, alt);
        int altDepth = low ? (int) (depth * vaf) : (int) (depth * (1 - vaf));
        evidence.ReadDepth = depth;
        evidence.RefSupport = depth - altDepth;
        evidence.AltSupport = altDepth;
        EvidenceByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(evidence);
    }

    private void createAndAddAmberSite(final HumanChromosome chromosome, final int position,
            String ref, String alt, final double gnomadFrequencyOfSite)
    {
        AmberSite site = new AmberSite(V38.versionedChromosome(chromosome), position, ref, alt, false, gnomadFrequencyOfSite);
        SitesByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(site);
    }

    private static int altDepthForGermlinePoint(double gnomadFrequency, int depth)
    {
        int altCount = 0;
        for(int i = 0; i < 2; i++)
        {
            if(Math.random() < gnomadFrequency)
            {
                altCount++;
            }
        }
        if(altCount == 0)
        {
            return 0;
        }
        if(altCount == 2)
        {
            return depth;
        }
        return depth / 2;
    }

    private static Pair<String, String> randomRefAndAltBases()
    {
        List<String> bases = new ArrayList<>(List.of("A", "C", "G", "T"));
        Collections.shuffle(bases);
        return Pair.of(bases.get(0), bases.get(1));
    }

    private PositionEvidence pe(String ref, String alt, double vaf)
    {
        final PositionEvidence pe = new PositionEvidence("1", 1000, ref, alt);
        pe.ReadDepth = 1000;
        pe.AltSupport = (int) (1000 * vaf);
        pe.RefSupport = 1000 - pe.AltSupport;
        return pe;
    }
}
