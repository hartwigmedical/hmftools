package com.hartwig.hmftools.amber.e2e;

import static com.hartwig.hmftools.amber.AmberConfig.LOCI_FILE;
import static com.hartwig.hmftools.amber.AmberConfig.USE_OLD_SEGMENTER;
import static com.hartwig.hmftools.amber.AmberConstants.TARGET_REGION_SITE_BUFFER;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.BLACKLISTED_SITES;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.amber.AmberApplication;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionImpl;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class RunAmberTest
{
    private String TumorSample;
    private String ReferenceSample;
    private File TumorBamFile;
    private File ReferenceBamFile;
    private File OutputDir;
    private Multimap<Chromosome, AmberBAF> Results;

    @Before
    public void setup() throws IOException
    {
        TumorSample = null;
        ReferenceSample = null;
        TumorBamFile = null;
        ReferenceBamFile = null;
        var TempDir = new File("/Users/timlavers/work/junk/basura"); //Files.createTempDirectory("amber").toFile();
        OutputDir = new File(TempDir, "output");
        //noinspection ResultOfMethodCallIgnored
        OutputDir.mkdirs();
        FileUtils.cleanDirectory(OutputDir);
        Results = null;
    }

    @Test
    public void twoChromosomes() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmberInWholeGenomeMode(scenario, false);
        // Check that baf.tsv file.
        scenario.checkResults(Results);

        // Check the segmentation file.
        String segmentsFile = PCFFile.generateBAFFilename(OutputDir.getAbsolutePath(), TumorSample);
        ListMultimap<Chromosome, PCFPosition> pcfData = PCFFile.readPositions(WINDOW_SIZE, PCFSource.TUMOR_BAF, segmentsFile);
        assertEquals(2, pcfData.keySet().size());
        List<PCFPosition> chr1Positions = pcfData.get(_1);
        // The R program that does segmentation puts a couple of spurious segments
        // at the start of each chromosome.
        assertEquals(9, chr1Positions.size());
        assertEquals(49001, chr1Positions.get(0).Position);
    }

    @Test
    public void newSegmenter() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmberInWholeGenomeMode(scenario, true);

        // Check the segmentation file.
        String segmentsFile = PCFFile.generateBAFFilename(OutputDir.getAbsolutePath(), TumorSample);
        ListMultimap<Chromosome, PcfSegment> pcfData = PCFFile.readPcfFile(segmentsFile);
        assertEquals(2, pcfData.keySet().size());
        List<PcfSegment> chr1Positions = pcfData.get(_1);
        assertEquals(3, chr1Positions.size());
        assertEquals(49554, chr1Positions.get(0).start());
    }

    @Test
    public void targetedMode() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmberInTargetedMode(scenario, "PanelA.bed", null);
        assertEquals(2, Results.keySet().size());
        ListMultimap<Chromosome, GenomeRegion> targetRegionsMap = ArrayListMultimap.create();
        targetRegionsMap.put(_1, new GenomeRegionImpl(_1, V38, 52000 - TARGET_REGION_SITE_BUFFER, 60000 + TARGET_REGION_SITE_BUFFER));
        targetRegionsMap.put(_2, new GenomeRegionImpl(_2, V38, 11000 - TARGET_REGION_SITE_BUFFER, 12000 + TARGET_REGION_SITE_BUFFER));
        scenario.checkResultsForTargetedMode(Results, targetRegionsMap);
    }

    @Test
    public void blacklistedPointsInTargetedMode() throws Exception
    {
        AmberScenario scenario = new AmberScenario("TwoChromosomes");
        runAmberInTargetedMode(scenario, "PanelA.bed", "BlacklistA.tsv");
        assertEquals(2, Results.keySet().size());
        ListMultimap<Chromosome, GenomeRegion> targetRegionsMap = ArrayListMultimap.create();
        targetRegionsMap.put(_1, new GenomeRegionImpl(_1, V38, 52000 - TARGET_REGION_SITE_BUFFER, 60000 + TARGET_REGION_SITE_BUFFER));
        targetRegionsMap.put(_2, new GenomeRegionImpl(_2, V38, 11000 - TARGET_REGION_SITE_BUFFER, 12000 + TARGET_REGION_SITE_BUFFER));
        ListMultimap<Chromosome, GenomePosition> excludedPositionsMap = ArrayListMultimap.create();
        excludedPositionsMap.put(_1, new GenomePositionImpl(_1, V38, 56381));
        excludedPositionsMap.put(_1, new GenomePositionImpl(_1, V38, 59040));
        excludedPositionsMap.put(_2, new GenomePositionImpl(_2, V38, 11486));
        excludedPositionsMap.put(_2, new GenomePositionImpl(_2, V38, 11607));
        excludedPositionsMap.put(_2, new GenomePositionImpl(_2, V38, 11834));

        scenario.checkResultsForTargetedModeWithExcludedPositions(Results, targetRegionsMap, excludedPositionsMap);
    }

    @Test
    public void handleEmptyData() throws Exception
    {
        AmberScenario scenario = new AmberScenario("NoBafs");
        runAmberInWholeGenomeMode(scenario, true);

        // Check the segmentation file.
        String segmentsFile = PCFFile.generateBAFFilename(OutputDir.getAbsolutePath(), TumorSample);
        ListMultimap<Chromosome, PcfSegment> pcfData = PCFFile.readPcfFile(segmentsFile);
        assertEquals(0, pcfData.keySet().size());
    }

    private void runAmberInTargetedMode(AmberScenario scenario, String panelFileName, String blacklistFileName) throws Exception
    {
        runAmber(scenario, false, panelFileName, blacklistFileName);
    }

    private void runAmberInWholeGenomeMode(AmberScenario scenario, boolean useNewSegmenter) throws Exception
    {
        runAmber(scenario, useNewSegmenter, null, null);
    }

    private void runAmber(AmberScenario scenario, boolean useNewSegmenter, String panelFileName, String blacklistFileName) throws Exception
    {
        File sitesFile = scenario.createAmberLocationsFile(OutputDir);
        TumorBamFile = scenario.getTumorBamFile();
        TumorSample = scenario.getTumorSampleName();
        File panelFile = panelFileName == null ? null : scenario.getTestDataFile(panelFileName);
        File blacklistFile = blacklistFileName == null ? null : scenario.getTestDataFile(blacklistFileName);

        int argCount = 6;
        if(TumorBamFile != null)
        {
            argCount += 4;
        }
        if(ReferenceBamFile != null)
        {
            argCount += 4;
        }
        if(panelFile != null)
        {
            Preconditions.checkArgument(panelFile.exists(), "Panel file does not exist: " + panelFile.getAbsolutePath());
            Preconditions.checkArgument(panelFile.isFile(), "Panel file is not a file: " + panelFile.getAbsolutePath());
            argCount += 2;
        }
        if(blacklistFileName != null)
        {
            Preconditions.checkNotNull(panelFile, "Blacklist file specified without panel file");
            Preconditions.checkArgument(blacklistFile.exists(), "Blacklist file does not exist: " + blacklistFile.getAbsolutePath());
            Preconditions.checkArgument(blacklistFile.isFile(), "Blacklist file is not a file: " + blacklistFile.getAbsolutePath());
            argCount += 2;
        }
        if(!useNewSegmenter)
        {
            argCount += 1;
        }
        String[] args = new String[argCount];
        int index = 0;
        args[index++] = String.format("-%s", LOCI_FILE);
        args[index++] = String.format("%s", sitesFile.getAbsolutePath());
        args[index++] = String.format("-%s", OUTPUT_DIR);
        args[index++] = String.format("%s", OutputDir.getAbsolutePath());
        args[index++] = String.format("-%s", REF_GENOME_VERSION);
        args[index++] = String.format("%s", "38");
        if(TumorBamFile != null)
        {
            args[index++] = String.format("-%s", TUMOR);
            args[index++] = String.format("%s", TumorSample);
            args[index++] = String.format("-%s", TUMOR_BAM);
            args[index++] = String.format("%s", TumorBamFile.getAbsolutePath());
        }
        if(ReferenceBamFile != null)
        {
            args[index++] = String.format("-%s", REFERENCE);
            args[index++] = String.format("%s", ReferenceSample);
            args[index++] = String.format("-%s", REFERENCE_BAM);
            args[index++] = String.format("%s", ReferenceBamFile.getAbsolutePath());
        }
        if(panelFile != null)
        {
            args[index++] = String.format("-%s", TARGET_REGIONS_BED);
            args[index++] = String.format("%s", panelFile.getAbsolutePath());
        }
        if(blacklistFileName != null)
        {
            args[index++] = String.format("-%s", BLACKLISTED_SITES);
            args[index++] = String.format("%s", blacklistFile.getAbsolutePath());
        }
        if(!useNewSegmenter)
        {
            args[index] = String.format("-%s", USE_OLD_SEGMENTER);
        }

        AmberApplication.main(args);

        File bafFile = new File(OutputDir, TumorSample + ".amber.baf.tsv.gz");
        Results = AmberBAFFile.read(bafFile.getAbsolutePath(), TumorSample != null);
    }
}
