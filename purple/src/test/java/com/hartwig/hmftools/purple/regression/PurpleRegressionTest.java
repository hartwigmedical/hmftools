package com.hartwig.hmftools.purple.regression;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.purple.PurpleApplication;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@SuppressWarnings("DataFlowIssue")
@Ignore("Manual run only")
public class PurpleRegressionTest
{
    private File WorkingDir = new File("/Users/timlavers/work/batches/2025/12/19/1");
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
    public void colo829Test() throws Exception
    {
        tumor = "COLO829v003T";
        reference = "COLO829v003R";
        File inputs = new File("/Users/timlavers/work/scratch/" + tumor + ".inputs.zip");
        Unzipper.unzipInto(inputs, InputDir);
        runPurple();
        File configuredOutputs = new File("/Users/timlavers/work/scratch/" + tumor + ".outputs.zip");
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

    private File makeWorkingSubDirectory(String name) throws IOException
    {
        File subDir = new File(WorkingDir, name);
        subDir.mkdirs();
        FileUtils.cleanDirectory(subDir);
        return subDir;
    }

    private void checkResults() throws IOException
    {
        File[] configuredResults = ConfiguredResultsDir.listFiles(File::isFile);
        File[] results = OutputDir.listFiles(File::isFile);
        Assert.assertEquals(configuredResults.length, results.length);
        boolean result = true;
        for (File configuredResult : configuredResults)
        {
            File resultFile = new File(OutputDir, configuredResult.getName());
            Assert.assertTrue(resultFile.exists());
            Assert.assertTrue(resultFile.isFile());
            if (resultFile.getName().contains(".gz"))
            {
                // Compare gz files by unzipping them and comparing the results line-by-line.
                result &= compareGzFiles(resultFile, configuredResult);
                continue;
            }
            if (resultFile.getName().contains("version"))
            {
                // Skip the version file as it contains a timestamp.
                continue;
            }
            // Compare other files by finding the first differing byte.
            long mismatchIndex = Files.mismatch(configuredResult.toPath(), resultFile.toPath());
            if (mismatchIndex == -1)
            {
                System.out.println("Match: " + configuredResult.getName());
            }
            else
            {
                System.out.println("Mismatch at position: " + mismatchIndex + " in " + configuredResult.getName());
                result = false;
            }
        }
        Assert.assertTrue(result);
    }

    private boolean compareGzFiles(File resultFile, File configuredFile) throws IOException
    {
        List<String> resultLines = IOUtils.readLines(FileWriterUtils.createBufferedReader(resultFile.getAbsolutePath()));
        List<String> configuredLines = IOUtils.readLines(FileWriterUtils.createBufferedReader(resultFile.getAbsolutePath()));
        if (resultLines.size() != configuredLines.size())
        {
            System.out.println("Mismatch in number of lines in file: " + resultFile.getName());
            return false;
        }
        for (int i = 0; i < resultLines.size(); i++)
        {
            if (!resultLines.get(i).equals(configuredLines.get(i)))
            {
                System.out.println("Mismatch at line " + i + " in file: " + resultFile.getName());
                return false;
            }
        }
        return true;
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
                "-somatic_sv_vcf", somatic_sv_vcf,
                "-germline_sv_vcf", germline_sv_vcf,
                "-somatic_vcf", somatic_vcf,
                "-germline_vcf", germline_vcf,
                "-driver_gene_panel", driver_gene_panel,
                "-amber", amber,
                "-cobalt", cobalt,
                "-ref_genome", ref_genome,
                "-ref_genome_version", ref_genome_version,
                "-gc_profile", gc_profile,
                "-somatic_hotspots", somatic_hotspots,
                "-germline_hotspots", germline_hotspots,
                "-ensembl_data_dir", ensembl_data_dir,
                "-germline_del_freq_file", germline_del_freq_file,
                "-threads", "16",
                //                "-circos", "/Users/timlavers/work/othercode/circos-0.69-9/bin/circos",
                "-output_dir", OutputDir.getAbsolutePath(), };
        PurpleApplication.main(args);
    }
}
