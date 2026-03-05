package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.COORDS_38;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
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

import org.junit.Before;
import org.junit.Test;

public class TumorOnlyPurityAnalysisTest extends PurityTestBase
{
    private static final String TUMOR_ID = "tumor123";
    private Map<Chromosome, SortedSet<AmberSite>> sitesByChromosome;
    private Map<Chromosome, SortedSet<PositionEvidence>> evidenceByChromosome;
    private File outputDir;
    private TumorOnlyPurityAnalysis analysis;

    @Before
    public void setup() throws Exception
    {
        sitesByChromosome = new HashMap<>();
        evidenceByChromosome = new HashMap<>();
        outputDir = Files.createTempDirectory("amber").toFile();
        outputDir.deleteOnExit();
    }

    @Test
    public void noPeaks() throws Exception
    {
        createBaselineReadings();
        createAnalysis();
        assertEquals(0.04, analysis.cutoff(), 0.001);
        assertTrue(analysis.contaminationPeaks().isEmpty());
    }

    @Test
    public void onePeakAboveMinimumCutoff() throws Exception
    {
        createBaselineReadings();
        final double vaf = 0.21;
        createContaminationPoints(vaf);
        createAnalysis();
        assertEquals(TumorOnlyPurityAnalysis.MIN_CUTOFF, analysis.cutoff(), 0.001);
        assertEquals(0, analysis.contaminationPeaks().size());
        assertEquals(1, analysis.copyNumberPeaks().size());
        assertEquals(vaf, analysis.copyNumberPeaks().get(0).vaf(), 0.01);
    }

    @Test
    public void onePeakBelowMinimumCutoff() throws Exception
    {
        createBaselineReadings();
        final double vaf = 0.1;
        // Add a 20 point copy number peak on 1P.
        createContaminationPoints(vaf);
        createAnalysis();
        assertEquals(vaf / 3, analysis.cutoff(), 0.01);
        assertEquals(0, analysis.contaminationPeaks().size());
        //        assertEquals(1, analysis.copyNumberPeaks().size());
        assertEquals(vaf, analysis.copyNumberPeaks().get(0).vaf(), 0.01);
    }

    private void createContaminationPoints(final double vaf)
    {
        // Add a 20 point copy number peak on 1P.
        int position = startPosition(new ChrArm(HumanChromosome._1, Arm.P)) + 50; // offset so as not to interfere with baseline points
        boolean low = false;
        for(int i = 0; i < 20; i++)
        {
            addPeakPoint(HumanChromosome._1, position, 1000, vaf, 0.5, low);
            position += 1000;
            low = !low;
        }
    }

    private void createBaselineReadings()
    {
        for(ChrArm chrArm : chromosomeArms())
        {
            int position = startPosition(chrArm);
            for(int i = 0; i < 100; i++)
            {
                addBaselinePoint(chrArm.chromosome(), position, 1000, 0.5);
                position += 1000;
            }
        }
    }

    private void createAnalysis() throws Exception
    {
        List<PositionEvidence> evidence = new ArrayList<>();
        for(Chromosome chromosome : evidenceByChromosome.keySet())
        {
            evidence.addAll(evidenceByChromosome.get(chromosome));
        }
        ListMultimap<Chromosome, AmberSite> sites = ArrayListMultimap.create();
        for(Chromosome chromosome : sitesByChromosome.keySet())
        {
            sites.putAll(chromosome, sitesByChromosome.get(chromosome));
        }
        PurityAnalysisConfig config = new PurityAnalysisConfig(TUMOR_ID, V38, outputDir.getAbsolutePath(), 8);
        analysis = new TumorOnlyPurityAnalysis(evidence, sites, config);
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
        createAndAddAmberSite(chromosome, position, gnomadFrequencyOfSite);

        PositionEvidence evidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, "A", "G");
        int altDepth = altDepthForNonContaminationPoint(gnomadFrequencyOfSite, depth);
        evidence.ReadDepth = depth;
        evidence.RefSupport = depth - altDepth;
        evidence.AltSupport = altDepth;
        evidenceByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(evidence);
    }

    private void addPeakPoint(HumanChromosome chromosome, int position, int depth, double vaf, double gnomadFrequencyOfSite, boolean low)
    {
        createAndAddAmberSite(chromosome, position, gnomadFrequencyOfSite);

        PositionEvidence evidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, "A", "G");
        int altDepth = low ? (int) (depth * vaf) : (int) (depth * (1 - vaf));
        evidence.ReadDepth = depth;
        evidence.RefSupport = depth - altDepth;
        evidence.AltSupport = altDepth;
        evidenceByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(evidence);
    }

    private void createAndAddAmberSite(final HumanChromosome chromosome, final int position, final double gnomadFrequencyOfSite)
    {
        AmberSite site = new AmberSite(V38.versionedChromosome(chromosome), position, "A", "G", false, gnomadFrequencyOfSite);
        sitesByChromosome.computeIfAbsent(chromosome, k -> new TreeSet<>()).add(site);
    }

    private static int altDepthForNonContaminationPoint(double gnomadFrequency, int depth)
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
}
