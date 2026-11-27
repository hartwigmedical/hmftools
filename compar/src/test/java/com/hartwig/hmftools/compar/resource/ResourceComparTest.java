package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.compar.resource.ResourceComparConfig.NEW_RESOURCES_DIR;
import static com.hartwig.hmftools.compar.resource.ResourceComparConfig.OLD_RESOURCES_DIR;
import static com.hartwig.hmftools.compar.resource.ResourceTestUtils.tsv;
import static com.hartwig.hmftools.compar.resource.ResourceTestUtils.writeFile;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ResourceComparTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    private File mOldDir;
    private File mNewDir;
    private File mOutputDir;

    private static final String OUTPUT_FILE = "compar_resources.tsv";

    private static final String HEADER = tsv("Category", "File", "MismatchType", "Key", "Differences", "OldValues", "NewValues");

    @Before
    public void setUp() throws IOException
    {
        mOldDir = temp.newFolder("old");
        mNewDir = temp.newFolder("new");
        mOutputDir = temp.newFolder("output");
    }

    @Test
    public void pairsFilesByNameAcrossFolders() throws IOException
    {
        writeBoth("DriverGenePanel.38.tsv", tsv("gene", "reportMissense"), tsv("KRAS", "true"));
        writeFile(mOldDir, "DriverGenePanel.38.tsv.dvc", "md5: 1");
        writeFile(mOldDir, "README.md", "notes");
        writeFile(mOldDir, "driver_genes.msk.37.tsv", tsv("gene", "reportMissense"), tsv("ALK", "true"));
        writeFile(mOldDir, "panel/tso500/driver_genes.tso500.38.tsv", tsv("gene", "reportMissense"), tsv("TP53", "true"));
        writeFile(mNewDir, "panel/tso500-v2/driver_genes.tso500.38.tsv", tsv("gene", "reportMissense"), tsv("TP53", "true"));
        writeFile(mNewDir, "known_fusions.38.bedpe", tsv("1", "100", "101", "2", "200", "201", "BCR_ABL1", "1", "+", "-"));

        assertTrue(runResourceCompar());

        assertEquals(List.of(
                HEADER,
                tsv("DRIVER_GENE_PANEL", "driver_genes.msk.37.tsv", "OLD_ONLY", "", "", "", ""),
                tsv("KNOWN_FUSION_BEDPE", "known_fusions.38.bedpe", "NEW_ONLY", "", "", "", "")),
                readOutput());

        assertEquals(List.of(OUTPUT_FILE), outputFiles());
    }

    @Test
    public void comparesRowsByKey() throws IOException
    {
        String header = tsv("gene", "reportMissense", "amplificationRatio", "reportGermlineVariant", "reportPGX");

        writeFile(mOldDir, "DriverGenePanel.38.tsv", header,
                tsv("KRAS", "true", "3", "NONE", "false"),
                tsv("ALK", "true", "3", "NONE", "false"),
                tsv("EGFR", "true", "3", "NONE", "false"),
                tsv("BRAF", "true", "3", "ANY", "false"),
                tsv("BRAF", "false", "3", "ANY", "false"),
                tsv("TP53", "true", "3", "false", "false"),
                tsv("PTEN", "true", "3", "true", "false"));

        writeFile(mNewDir, "DriverGenePanel.38.tsv", header,
                tsv("KRAS", "true", "2.5", "NONE", "false"),
                tsv("EGFR", "true", "3", "NONE", "false"),
                tsv("EGFR", "true", "2", "NONE", "false"),
                tsv("BRAF", "true", "3", "ANY", "false"),
                tsv("BRAF", "false", "3", "ANY", "false"),
                tsv("TP53", "TRUE", "3", "NONE", "true"),
                tsv("PTEN", "true", "3", "ANY", "false"),
                tsv("MYC", "true", "3", "ANY", "false"));

        assertTrue(runResourceCompar());

        String file = "DriverGenePanel.38.tsv";

        assertEquals(List.of(
                HEADER,
                tsv("DRIVER_GENE_PANEL", file, "VALUE", "KRAS", "amplificationRatio(3.0/2.5)",
                        panelValues("3.0", "NONE", "false"), panelValues("2.5", "NONE", "false")),
                tsv("DRIVER_GENE_PANEL", file, "OLD_ONLY", "ALK", "", panelValues("3.0", "NONE", "false"), ""),
                tsv("DRIVER_GENE_PANEL", file, "DUPLICATE_KEY", "EGFR", "rows(1/2)", panelValues("3.0", "NONE", "false"), ""),
                tsv("DRIVER_GENE_PANEL", file, "DUPLICATE_KEY", "EGFR", "rows(1/2)", "", panelValues("3.0", "NONE", "false")),
                tsv("DRIVER_GENE_PANEL", file, "DUPLICATE_KEY", "EGFR", "rows(1/2)", "", panelValues("2.0", "NONE", "false")),
                tsv("DRIVER_GENE_PANEL", file, "VALUE", "TP53", "reportPGX(false/true)",
                        panelValues("3.0", "NONE", "false"), panelValues("3.0", "NONE", "true")),
                tsv("DRIVER_GENE_PANEL", file, "NEW_ONLY", "MYC", "", "", panelValues("3.0", "ANY", "false"))),
                readOutput());
    }

    @Test
    public void comparesMissingPanelColumnsAsDefaults() throws IOException
    {
        writeFile(mOldDir, "DriverGenePanel.38.tsv",
                tsv("gene", "reportLOH", "amplificationRatio"), tsv("BRCA1", "true", "3"), tsv("KRAS", "false", "3"));

        writeFile(mNewDir, "DriverGenePanel.38.tsv",
                tsv("gene", "reportLoh", "amplificationRatio"), tsv("BRCA1", "true", "3.0"), tsv("KRAS", "false", "3"));

        assertTrue(runResourceCompar());

        String file = "DriverGenePanel.38.tsv";

        String defaultValues = "reportHetDeletion=false;reportLoh=%s;hetDeletionThreshold=0.6;amplificationRatio=3.0;"
                + "reportGermlineAmplification=false;reportHighExpression=false;reportLowExpression=false;reportNovelSpliceJunction=false";

        assertEquals(List.of(
                HEADER,
                tsv("DRIVER_GENE_PANEL", file, "VALUE", "", "columns(reportLOH/reportLoh)", "", ""),
                tsv("DRIVER_GENE_PANEL", file, "VALUE", "BRCA1", "reportLoh(false/true)",
                        String.format(defaultValues, "false"), String.format(defaultValues, "true"))),
                readOutput());
    }

    @Test
    public void comparesColumnSetsAtFileLevel() throws IOException
    {
        writeFile(mOldDir, "driver_genes.msk.37.tsv", tsv("gene", "reportMissense", "reportPGX"), tsv("ALK", "true", "false"));
        writeFile(mNewDir, "driver_genes.msk.37.tsv", tsv("gene", "reportMissense", "reportLoh"), tsv("ALK", "true", "false"));

        writeFile(mOldDir, "driver_genes.tso500.38.tsv", tsv("gene", "reportMissense", "amplificationRatio"), tsv("KRAS", "true", "3"));
        writeFile(mNewDir, "driver_genes.tso500.38.tsv", tsv("gene", "amplificationRatio", "reportMissense"), tsv("KRAS", "3", "true"));

        assertTrue(runResourceCompar());

        assertEquals(List.of(
                HEADER,
                tsv("DRIVER_GENE_PANEL", "driver_genes.msk.37.tsv", "VALUE", "", "columns(reportPGX/reportLoh)", "", "")),
                readOutput());
    }

    @Test
    public void comparesInfoFieldsByKey() throws IOException
    {
        String vcfHeader = tsv("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO");

        writeFile(mOldDir, "clinvar.38.vcf.gz", "##fileformat=VCFv4.1", vcfHeader,
                tsv("chr1", "100", "1", "A", "G", ".", ".", "CLNSIG=Benign;GENEINFO=GENE1:1"),
                tsv("chr1", "200", "2", "C", "T", ".", ".", "CLNSIG=Pathogenic;GENEINFO=GENE2:2"),
                tsv("chr1", "300", "3", "G", "A", ".", ".", "."));

        writeFile(mNewDir, "clinvar.38.vcf.gz", "##fileformat=VCFv4.1", vcfHeader,
                tsv("1", "100", "1", "A", "G", ".", ".", "CLNSIG=Benign;GENEINFO=GENE1:9"),
                tsv("1", "200", "2", "C", "T", ".", ".", "CLNSIG=Likely_pathogenic;CLNSIGCONF=x;GENEINFO=GENE2:2"),
                tsv("1", "300", "3", "G", "A", ".", ".", "."),
                tsv("1", "400", "4", "T", "C", ".", ".", "CLNSIG=Benign;CLNVC"));

        String hotspotHeader = tsv("Chromosome", "Position", "Ref", "Alt", "Info");

        writeFile(mOldDir, "KnownHotspots.somatic.38.tsv.gz", hotspotHeader,
                tsv("chr1", "500", "T", "C", "input=a"),
                tsv("chr2", "700", "G", "T", ""));

        writeFile(mNewDir, "KnownHotspots.somatic.38.tsv.gz", hotspotHeader,
                tsv("chr1", "500", "T", "C", "input=b"));

        assertTrue(runResourceCompar());

        String clinvar = "clinvar.38.vcf.gz";
        String hotspots = "KnownHotspots.somatic.38.tsv.gz";

        assertEquals(List.of(
                HEADER,
                tsv("SOMATIC_HOTSPOT_TSV", hotspots, "OLD_ONLY", "chr2;700;G;T", "", "input=", ""),
                tsv("CLINVAR", clinvar, "VALUE", "", "columns(/CLNSIGCONF;CLNVC)", "", ""),
                tsv("CLINVAR", clinvar, "VALUE", "1;200;C;T", "CLNSIG(Pathogenic/Likely_pathogenic)",
                        "CLNSIG=Pathogenic;GENEINFO=GENE2:2", "CLNSIG=Likely_pathogenic;CLNSIGCONF=x;CLNVC=;GENEINFO=GENE2:2"),
                tsv("CLINVAR", clinvar, "NEW_ONLY", "1;400;T;C", "", "", "CLNSIG=Benign;CLNSIGCONF=;CLNVC=true;GENEINFO=")),
                readOutput());
    }

    @Test
    public void comparesCsvAndHeaderlessFiles() throws IOException
    {
        String csvHeader = "Type,FiveGene,ThreeGene,CancerTypes,PubMedId,KnownExonTranscript,KnownExonUpRange,KnownExonDownRange,"
                + "HighImpactPromiscuous,Overrides";

        writeFile(mOldDir, "known_fusion_data.38.csv", csvHeader,
                "KNOWN_PAIR,BCR,ABL1,Leukaemia,1,,,,FALSE,",
                "PROMISCUOUS_3,,RET,Thyroid,2,ENST1,,11;12,false,");

        writeFile(mNewDir, "known_fusion_data.38.csv", csvHeader,
                "KNOWN_PAIR,BCR,ABL1,Leukaemia,1,,,,false,",
                "PROMISCUOUS_3,,RET,Thyroid,2,ENST1,,11;13,false,");

        writeFile(mOldDir, "known_fusions.38.bedpe",
                tsv("1", "100", "101", "2", "200", "201", "BCR_ABL1", "1", "+", "-"),
                tsv("3", "300", "301", "4", "400", "401", "EML4_ALK", "1", "+", "+"));

        writeFile(mNewDir, "known_fusions.38.bedpe",
                tsv("1", "100", "101", "2", "200", "201", "BCR_ABL1", "5", "+", "-"),
                tsv("3", "300", "301", "4", "400", "401", "EML4_ALK", "1", "-", "+"));

        assertTrue(runResourceCompar());

        String csv = "known_fusion_data.38.csv";
        String bedpe = "known_fusions.38.bedpe";

        String fusionValues = "CancerTypes=Thyroid;PubMedId=2;KnownExonTranscript=ENST1;KnownExonUpRange=;KnownExonDownRange=%s;"
                + "HighImpactPromiscuous=false;Overrides=";

        assertEquals(List.of(
                HEADER,
                tsv("KNOWN_FUSION_DATA", csv, "VALUE", "PROMISCUOUS_3;;RET", "KnownExonDownRange(11;12/11;13)",
                        String.format(fusionValues, "11;12"), String.format(fusionValues, "11;13")),
                tsv("KNOWN_FUSION_BEDPE", bedpe, "OLD_ONLY", "3;300;301;4;400;401;+;+", "", "name=EML4_ALK;score=1", ""),
                tsv("KNOWN_FUSION_BEDPE", bedpe, "NEW_ONLY", "3;300;301;4;400;401;-;+", "", "", "name=EML4_ALK;score=1")),
                readOutput());
    }

    @Test
    public void comparesPonFiles() throws IOException
    {
        String header = tsv("Chromosome", "Position", "Ref", "Alt", "SampleCount", "MaxSampleReadCount", "TotalReadCount");

        writeFile(mOldDir, "hmf_wgs_sage_pon_1000.38.tsv.gz", header,
                tsv("chr1", "100", "A", "G", "4", "5", "18"),
                tsv("chr3", "99999997", "GAAAAGA", "G", "3", "2", "6"),
                tsv("chr3", "99999997", "GAAAAGA", "G", "3", "2", "6"));

        writeFile(mNewDir, "hmf_wgs_sage_pon_1000.38.tsv.gz", header,
                tsv("chr1", "100", "A", "G", "5", "5", "19"),
                tsv("chr3", "99999997", "GAAAAGA", "G", "3", "2", "6"));

        writeFile(mOldDir, "sv_pon.38.bedpe.gz", tsv("1", "100", "101", "1", "500", "501", ".", "2", "+", "-"));
        writeFile(mNewDir, "sv_pon.38.bedpe.gz", tsv("1", "100", "101", "1", "500", "501", ".", "3", "+", "-"));

        assertTrue(runResourceCompar());

        String sagePon = "hmf_wgs_sage_pon_1000.38.tsv.gz";
        String svPon = "sv_pon.38.bedpe.gz";
        String duplicateValues = "SampleCount=3;MaxSampleReadCount=2;TotalReadCount=6";

        assertEquals(List.of(
                HEADER,
                tsv("SAGE_PON", sagePon, "VALUE", "chr1;100;A;G", "SampleCount(4/5);TotalReadCount(18/19)",
                        "SampleCount=4;MaxSampleReadCount=5;TotalReadCount=18", "SampleCount=5;MaxSampleReadCount=5;TotalReadCount=19"),
                tsv("SAGE_PON", sagePon, "DUPLICATE_KEY", "chr3;99999997;GAAAAGA;G", "rows(2/1)", duplicateValues, ""),
                tsv("SAGE_PON", sagePon, "DUPLICATE_KEY", "chr3;99999997;GAAAAGA;G", "rows(2/1)", duplicateValues, ""),
                tsv("SAGE_PON", sagePon, "DUPLICATE_KEY", "chr3;99999997;GAAAAGA;G", "rows(2/1)", "", duplicateValues),
                tsv("SV_PON", svPon, "VALUE", "1;100;101;1;500;501;+;-", "score(2/3)", "score=2;name=.", "score=3;name=.")),
                readOutput());
    }

    @Test
    public void reportsInvalidFiles() throws IOException
    {
        writeFile(mOldDir, "DriverGenePanel.38.tsv", tsv("gene", "reportMissense"), tsv("KRAS", "true"), tsv("ALK", "true"));
        writeFile(mNewDir, "DriverGenePanel.38.tsv", tsv("name", "reportMissense"), tsv("KRAS", "false"));

        writeFile(mOldDir, "known_fusion_data.38.csv", "Type,FiveGene,ThreeGene", "KNOWN_PAIR,BCR");
        writeFile(mNewDir, "known_fusion_data.38.csv", "Type,FiveGene,ThreeGene", "KNOWN_PAIR,BCR,ABL1");

        Path oldVcf = writeFile(mOldDir, "clinvar.38.vcf.gz", "##fileformat=VCFv4.1",
                tsv("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
                tsv("1", "100", "1", "A", "G", ".", ".", "CLNSIG=Benign")).toPath();
        byte[] vcfBytes = Files.readAllBytes(oldVcf);
        Files.write(oldVcf, Arrays.copyOf(vcfBytes, vcfBytes.length / 2));
        Files.write(mNewDir.toPath().resolve("clinvar.38.vcf.gz"), Arrays.copyOf(vcfBytes, vcfBytes.length / 3));

        writeFile(mOldDir, "hmf_wgs_sage_pon_1000.38.tsv.gz");
        Files.writeString(mNewDir.toPath().resolve("hmf_wgs_sage_pon_1000.38.tsv.gz"), "not gzip");

        assertTrue(runResourceCompar());

        assertEquals(List.of(
                HEADER,
                tsv("DRIVER_GENE_PANEL", "DriverGenePanel.38.tsv", "INVALID_NEW", "", "", "", ""),
                tsv("KNOWN_FUSION_DATA", "known_fusion_data.38.csv", "INVALID_OLD", "", "", "", ""),
                tsv("CLINVAR", "clinvar.38.vcf.gz", "INVALID_BOTH", "", "", "", ""),
                tsv("SAGE_PON", "hmf_wgs_sage_pon_1000.38.tsv.gz", "INVALID_BOTH", "", "", "", "")),
                readOutput());
    }

    @Test
    public void failsWhenFileNameIsFoundTwice() throws IOException
    {
        writeFile(mOldDir, "a/DriverGenePanel.38.tsv", tsv("gene", "reportMissense"), tsv("KRAS", "true"));
        writeFile(mOldDir, "b/DriverGenePanel.38.tsv", tsv("gene", "reportMissense"), tsv("KRAS", "true"));

        assertFalse(runResourceCompar());
    }

    @Test
    public void mainWritesOutputFileWithOutputId()
    {
        ResourceCompar.main(new String[] {
                "-old_resources_dir", mOldDir.getPath(),
                "-new_resources_dir", mNewDir.getPath(),
                "-output_dir", mOutputDir.getPath(),
                "-output_id", "release" });

        assertEquals(List.of("compar_resources.release.tsv"), outputFiles());
    }

    private boolean runResourceCompar()
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Compar");
        ResourceComparConfig.addConfig(configBuilder);
        configBuilder.setValue(OLD_RESOURCES_DIR, mOldDir.getPath());
        configBuilder.setValue(NEW_RESOURCES_DIR, mNewDir.getPath());
        configBuilder.setValue(OUTPUT_DIR, mOutputDir.getPath());
        return new ResourceCompar(configBuilder).run();
    }

    private void writeBoth(final String relativePath, final String... lines) throws IOException
    {
        writeFile(mOldDir, relativePath, lines);
        writeFile(mNewDir, relativePath, lines);
    }

    private static String panelValues(final String amplificationRatio, final String germlineVariant, final String reportPGX)
    {
        return String.format(
                "reportMissense=true;reportHetDeletion=false;reportLoh=false;hetDeletionThreshold=0.6;amplificationRatio=%s;"
                        + "reportGermlineVariant=%s;reportGermlineAmplification=false;reportPGX=%s;reportHighExpression=false;"
                        + "reportLowExpression=false;reportNovelSpliceJunction=false",
                amplificationRatio, germlineVariant, reportPGX);
    }

    private List<String> readOutput() throws IOException
    {
        return Files.readAllLines(mOutputDir.toPath().resolve(OUTPUT_FILE));
    }

    private List<String> outputFiles()
    {
        String[] fileNames = mOutputDir.list();
        Arrays.sort(fileNames);
        return Arrays.asList(fileNames);
    }
}
