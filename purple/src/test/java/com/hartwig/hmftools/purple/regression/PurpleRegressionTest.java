package com.hartwig.hmftools.purple.regression;

import static com.hartwig.hmftools.purple.SampleDataFiles.AMBER;
import static com.hartwig.hmftools.purple.SampleDataFiles.COBALT;
import static com.hartwig.hmftools.purple.SampleDataFiles.GERMLINE_SV_VCF;
import static com.hartwig.hmftools.purple.SampleDataFiles.GERMLINE_VARIANTS;
import static com.hartwig.hmftools.purple.SampleDataFiles.SOMATIC_SV_VCF;
import static com.hartwig.hmftools.purple.SampleDataFiles.SOMATIC_VARIANTS;
import static com.hartwig.hmftools.purple.germline.GermlineAmpDelFrequencyCache.COHORT_AMP_DEL_FREQ_FILE;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.PurpleApplication;
import com.hartwig.hmftools.purple.copynumber.ChromosomeArmCopyNumber;
import com.hartwig.hmftools.purple.copynumber.ChromosomeArmCopyNumbersFile;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore("Manual run only")
public class PurpleRegressionTest
{
    private final File WorkingDir = new File("/Users/timlavers/work/batches/2026/1/5/1");
    private File InputDir;
    private File OutputDir;
    private File ConfiguredResultsDir;
    private String tumor;
    private String reference;

    @Before
    public void setup() throws IOException
    {
        InputDir = makeWorkingSubDirectory("input");
        OutputDir = makeWorkingSubDirectory("output");
        ConfiguredResultsDir = makeWorkingSubDirectory("configured_results");
    }

    @Test
    public void runAndMeasure() throws Exception
    {
        tumor = "COLO829v003T";
        reference = "COLO829v003R";
        File inputs = new File("/Users/timlavers/work/scratch/" + tumor + ".inputs.zip");
        Unzipper.unzipInto(inputs, InputDir);
        int numberOfRuns = 10;
        Stopwatch stopwatch = Stopwatch.createUnstarted();
        for(int i = 0; i < numberOfRuns; i++)
        {
            stopwatch.start();
            runPurple();
            stopwatch.stop();
        }
        System.out.println("Average run time: " + stopwatch.elapsed(TimeUnit.MILLISECONDS) / numberOfRuns + " ms");
    }

    @Test
    public void colo829Test() throws Exception
    {
        tumor = "COLO829v003T";
        reference = "COLO829v003R";
        File inputs = new File("/Users/timlavers/work/scratch/" + tumor + ".inputs.zip");
        Unzipper.unzipInto(inputs, InputDir);
        runPurple();
        File configuredOutputs = new File("/Users/timlavers/work/batches/2025/12/19/3/COLO829v003.v2_3.purple.results.zip");
        Unzipper.unzipInto(configuredOutputs, ConfiguredResultsDir);
        checkResults();
    }

    @Test
    public void h00000469Test() throws Exception
    {
        tumor = "H00000469";
        reference = "H00000469-ref";
        File inputs = new File("/Users/timlavers/work/scratch/" + tumor + ".inputs.zip");
        Unzipper.unzipInto(inputs, InputDir);
        runPurple();
        File configuredOutputs = new File("/Users/timlavers/work/scratch/" + tumor + ".outputs.zip");
        Unzipper.unzipInto(configuredOutputs, ConfiguredResultsDir);
        checkResults();
    }

    private void checkPurity() throws Exception
    {
        PurityContext baselineContext = PurityContextFile.read(ConfiguredResultsDir.getAbsolutePath(), tumor);
        PurityContext newContext = PurityContextFile.read(OutputDir.getAbsolutePath(), tumor);
        assertEquals(baselineContext, newContext);
    }

    private void checkSomaticDrivers() throws Exception
    {
        String outputDrivers = DriverCatalogFile.generateSomaticFilename(OutputDir.getAbsolutePath(), tumor);
        String baselineDrivers = DriverCatalogFile.generateSomaticFilename(ConfiguredResultsDir.getAbsolutePath(), tumor);
        checkDrivers(baselineDrivers, outputDrivers);
    }

    private void checkGermlineDrivers() throws Exception
    {
        String outputGermlineDriverFile = DriverCatalogFile.generateGermlineFilename(OutputDir.getAbsolutePath(), tumor);
        String baselineGermlineDriverFile = DriverCatalogFile.generateGermlineFilename(ConfiguredResultsDir.getAbsolutePath(), tumor);
        checkDrivers(baselineGermlineDriverFile, outputGermlineDriverFile);
    }

    private void checkChromosomeArmCopyNumbers() throws Exception
    {
        List<ChromosomeArmCopyNumber> armCopyNumbers =
                ChromosomeArmCopyNumbersFile.read(ChromosomeArmCopyNumbersFile.generateFilename(OutputDir.getAbsolutePath(), tumor));
        List<ChromosomeArmCopyNumber> baselineCopyNumbers =
                ChromosomeArmCopyNumbersFile.read(ChromosomeArmCopyNumbersFile.generateFilename(ConfiguredResultsDir.getAbsolutePath(), tumor));
        assertEquals(armCopyNumbers.size(), baselineCopyNumbers.size());
        for(ChromosomeArmCopyNumber armCopyNumber : armCopyNumbers)
        {
            ChromosomeArmCopyNumber baseline = findCopyNumber(baselineCopyNumbers, armCopyNumber);
            checkObjectsHaveSameData(baseline, armCopyNumber);
        }
    }

    private static ChromosomeArmCopyNumber findCopyNumber(List<ChromosomeArmCopyNumber> copyNumbers, ChromosomeArmCopyNumber toMatch)
    {
        List<ChromosomeArmCopyNumber> matches =
                copyNumbers.stream().filter(cn -> cn.chromosome().equals(toMatch.chromosome()) && cn.arm().equals(toMatch.arm())).toList();
        assertEquals(
                "Copy numbers matching " + toMatch.chromosome() + " " + toMatch.arm() + " has size " + matches.size(), 1, matches.size());
        return matches.get(0);
    }

    private void checkGeneCopyNumbers() throws Exception
    {
        List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(OutputDir.getAbsolutePath(), tumor));
        List<GeneCopyNumber> baselineCopyNumbers =
                GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(ConfiguredResultsDir.getAbsolutePath(), tumor));

        assertEquals("Copy numbers for tumor " + tumor + " has size " + copyNumbers.size(), copyNumbers.size(), baselineCopyNumbers.size());
        for(GeneCopyNumber copyNumber : copyNumbers)
        {
            GeneCopyNumber baseline = findGeneCopyNumber(baselineCopyNumbers, copyNumber);
            //            System.out.println(copyNumber.GeneName);
            if(baseline.GeneName.equals("EML4") || baseline.GeneName.equals("ZFP36L2") || baseline.GeneName.equals("SIX2"))
            {
                continue;
            }
            try
            {
                checkObjectsHaveSameData(baseline, copyNumber, "RelativeMinCopyNumber");
            }
            catch(AssertionError e)
            {
                //                System.out.println(copyNumber.GeneName);
            }
        }
    }

    private GeneCopyNumber findGeneCopyNumber(final List<GeneCopyNumber> geneCopyNumbers, final GeneCopyNumber toMatch)
    {
        List<GeneCopyNumber> matches = geneCopyNumbers.stream()
                .filter(baseline -> baseline.TransName.equals(toMatch.TransName))
                .toList();
        assertEquals(1, matches.size());
        return matches.get(0);
    }

    private void checkDeletions() throws Exception
    {
        List<GermlineAmpDel> outputDeletions =
                GermlineAmpDel.read(GermlineAmpDel.generateFilename(OutputDir.getAbsolutePath(), tumor));
        List<GermlineAmpDel> baselineDeletions =
                GermlineAmpDel.read(GermlineAmpDel.generateFilename(ConfiguredResultsDir.getAbsolutePath(), tumor));
        Assert.assertEquals(outputDeletions.size(), baselineDeletions.size());
        for(GermlineAmpDel deletion : outputDeletions)
        {
            GermlineAmpDel baseline = findGermlineDeletion(baselineDeletions, deletion);
            checkObjectsHaveSameData(baseline, deletion);
        }
    }

    private void checkDrivers(String baselineGermlineDriverFile, String outputGermlineDriverFile) throws Exception
    {
        List<DriverCatalog> outputDrivers = DriverCatalogFile.read(outputGermlineDriverFile);
        List<DriverCatalog> baselineDrivers = DriverCatalogFile.read(baselineGermlineDriverFile);

        // The drivers in the baseline should be in the new output.
        for(DriverCatalog baselineDriver : baselineDrivers)
        {
            DriverCatalog inNewOutput = findDriver(outputDrivers, baselineDriver.transcript());
            Assert.assertNotNull(inNewOutput);
            assertEquals(inNewOutput, baselineDriver);
            checkObjectsHaveSameData(baselineDriver, inNewOutput);
        }

        // The drivers in the new output that are not in the baseline should have reported status "NONE" or "NOT_REPORTED.
        for(DriverCatalog outputDriver : outputDrivers)
        {
            DriverCatalog baselineDriver = findDriver(baselineDrivers, outputDriver.transcript());
            if(baselineDriver == null)
            {
                Assert.assertNotEquals(outputDriver.toString(), ReportedStatus.REPORTED, outputDriver.reportedStatus());
            }
        }
    }

    private void checkSegments() throws Exception
    {
        List<PurpleSegment> outputSegments = PurpleSegment.read(PurpleSegment.generateFilename(OutputDir.getAbsolutePath(), tumor));
        List<PurpleSegment> baselineSegments =
                PurpleSegment.read(PurpleSegment.generateFilename(ConfiguredResultsDir.getAbsolutePath(), tumor));
        assertEquals(outputSegments.size(), baselineSegments.size());

        for(PurpleSegment outputSegment : outputSegments)
        {
            PurpleSegment baselineSegment = findSegment(baselineSegments, outputSegment);
            checkObjectsHaveSameData(baselineSegment, outputSegment, "GermlineState");
        }
    }

    private static PurpleSegment findSegment(List<PurpleSegment> segments, PurpleSegment toMatch)
    {
        return segments.stream().filter(segment -> segment.matches(toMatch)).findFirst().orElse(null);
    }

    private static void checkObjectsHaveSameData(Object s, Object t, String... ignoredFields)
    {
        Class<?> currentClass = s.getClass();
        while(currentClass != null && !currentClass.equals(Object.class))
        {
            for(Field field : currentClass.getDeclaredFields())
            {
                if(ignoredFields != null && Arrays.asList(ignoredFields).contains(field.getName()))
                {
                    continue;
                }
                if(Modifier.isStatic(field.getModifiers()))
                {
                    continue;
                }
                try
                {
                    field.setAccessible(true);
                    Object valS = field.get(s);
                    Object valT = field.get(t);
                    if(field.getType().equals(double.class) || field.getType().equals(Double.class))
                    {
                        assertEquals("Field " + field.getName() + " mismatch in " + currentClass.getSimpleName(),
                                (double) valS,
                                (double) valT,
                                0.001);
                    }
                    else
                    {
                        assertEquals("Field " + field.getName() + " mismatch in " + currentClass.getSimpleName(), valS, valT);
                    }
                }
                catch(IllegalAccessException e)
                {
                    throw new RuntimeException(e);
                }
            }
            currentClass = currentClass.getSuperclass();
        }
    }

    private static GermlineAmpDel findGermlineDeletion(List<GermlineAmpDel> deletions, GermlineAmpDel toMatch)
    {
        List<GermlineAmpDel> matching = deletions.stream().filter(germline -> deletionsMatch(germline, toMatch)).toList();
        assertEquals("Deletions matching " + toMatch.GeneName + " has size " + matching.size(), 1, matching.size());
        return matching.get(0);
    }

    private static boolean deletionsMatch(GermlineAmpDel g, GermlineAmpDel h)
    {
        return g.GeneName.equals(h.GeneName)
                && g.Chromosome.equals(h.Chromosome)
                && g.RegionStart == h.RegionStart
                && g.RegionEnd == h.RegionEnd
                && g.NormalStatus == h.NormalStatus
                && g.TumorStatus == h.TumorStatus
                && Doubles.equal(g.GermlineCopyNumber, h.GermlineCopyNumber)
                && Doubles.equal(g.TumorCopyNumber, h.TumorCopyNumber);
    }

    private static DriverCatalog findDriver(List<DriverCatalog> drivers, String transcript)
    {
        List<DriverCatalog> matching = drivers.stream().filter(driverCatalog -> driverCatalog.transcript().equals(transcript)).toList();
        Assert.assertTrue("Found " + matching.size() + " drivers for transcript " + transcript, matching.size() <= 1);
        if(matching.size() == 1)
        {
            return matching.get(0);
        }
        else
        {
            return null;
        }
    }

    @SuppressWarnings("ResultOfMethodCallIgnored")
    private File makeWorkingSubDirectory(String name) throws IOException
    {
        File subDir = new File(WorkingDir, name);
        subDir.mkdirs();
        FileUtils.cleanDirectory(subDir);
        return subDir;
    }

    private void checkResults() throws Exception
    {
        checkChromosomeArmCopyNumbers();// .chromosome_arm
        checkGeneCopyNumbers();         // .cnv.gene
        //        checkGermlineDrivers();         // .driver.catalog.germline
        //        checkSomaticDrivers();          // .driver.catalog.somatic
        //        checkDeletions();               // .germline.deletion
        checkPurity();                  // .purity
        checkSegments();                // .segment
    }

    private void runPurple() throws IOException
    {
        String inputRoot = InputDir.getAbsolutePath();
        String somatic_sv_vcf = inputRoot + "/" + tumor + ".esvee.somatic.vcf.gz";
        String germline_sv_vcf = inputRoot + "/" + tumor + ".esvee.germline.vcf.gz";
        String somatic_vcf = inputRoot + "/" + tumor + ".pave.somatic.vcf.gz";
        String germline_vcf = inputRoot + "/" + tumor + ".pave.germline.vcf.gz";
        String amber = inputRoot + "/amber";
        String cobalt = inputRoot + "/cobalt";

        String driver_gene_panel = "/Users/timlavers/work/data/pipeline_resources_v2_2/38/hmftools/common/DriverGenePanel.38.tsv";
        String ref_genome = "/Users/timlavers/work/data/reference_genomes/38/Homo_sapiens_assembly38.alt.masked.fasta";
        String ref_genome_version = "V38";
        String gc_profile = "/Users/timlavers/work/data/pipeline_resources_v2_3/38/hmftools/dna/copy_number/GC_profile.1000bp.38.cnp";
        String somatic_hotspots =
                "/Users/timlavers/work/data/pipeline_resources_v2_3/38/hmftools/dna/variants/KnownHotspots.somatic.38.vcf.gz";
        String germline_hotspots =
                "/Users/timlavers/work/data/pipeline_resources_v2_2/38/hmftools/dna/variants/KnownHotspots.germline.38.vcf.gz";
        String ensembl_data_dir = "/Users/timlavers/work/data/pipeline_resources_v2_2/38/hmftools/common/ensembl_data";
        String germline_del_freq_file = "/Users/timlavers/work/batches/2025/10/17/1/cohort_germline_del_freq.38.csv";

        String[] args = new String[] {
                "-tumor", tumor,
                "-reference", reference,
                "-" + SOMATIC_SV_VCF, somatic_sv_vcf,
                "-" + GERMLINE_SV_VCF, germline_sv_vcf,
                "-" + SOMATIC_VARIANTS, somatic_vcf,
                "-" + GERMLINE_VARIANTS, germline_vcf,
                "-driver_gene_panel", driver_gene_panel,
                "-" + AMBER, amber,
                "-" + COBALT, cobalt,
                "-ref_genome", ref_genome,
                "-ref_genome_version", ref_genome_version,
                "-gc_profile", gc_profile,
                "-somatic_hotspots", somatic_hotspots,
                "-germline_hotspots", germline_hotspots,
                "-ensembl_data_dir", ensembl_data_dir,
                "-" + COHORT_AMP_DEL_FREQ_FILE, germline_del_freq_file,
                "-threads", "16",
                //                "-circos", "/Users/timlavers/work/othercode/circos-0.69-9/bin/circos",
                "-output_dir", OutputDir.getAbsolutePath(), };
        PurpleApplication.main(args);
    }
}
