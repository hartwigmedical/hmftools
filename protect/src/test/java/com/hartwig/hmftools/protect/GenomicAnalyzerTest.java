package com.hartwig.hmftools.protect;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.junit.Test;

public class GenomicAnalyzerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String BACHELOR_TSV = BASE_DIRECTORY + "/bachelor/sample.reportable_germline_variants.tsv";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_VIRAL_INSERTION_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.drivers.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";

    @Test
    public void canRunOnTestRun() throws IOException {
        GenomicAnalyzer analyzer =
                new GenomicAnalyzer(ProtectTestFactory.loadTestActionabilityAnalyzer(), new GermlineReportingModel(Maps.newHashMap()));

        PatientTumorLocation patientTumorLocation = null;

        assertNotNull(analyzer.run("sample",
                patientTumorLocation,
                PURPLE_PURITY_TSV,
                PURPLE_QC_FILE,
                PURPLE_DRIVER_CATALOG_TSV,
                PURPLE_SOMATIC_VARIANT_VCF,
                BACHELOR_TSV,
                LINX_FUSIONS_TSV,
                LINX_BREAKEND_TSV,
                LINX_VIRAL_INSERTION_TSV,
                LINX_DRIVERS_TSV,
                CHORD_PREDICTION_TXT));
    }
}
