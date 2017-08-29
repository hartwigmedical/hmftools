package com.hartwig.hmftools.patientreporter.report;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.HmfReporterDataLoader;

import org.junit.Test;

import net.sf.jasperreports.engine.JRDataSource;

public class GenePanelDataSourceTest {

    @Test
    public void canCreateGenePanelFor3Genes() throws IOException, HartwigException {
        final String slicerPath = Resources.getResource("bed").getPath() + File.separator + "hmf_gene_panel.tsv";
        final String drupFilterPath = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
        final String cosmicPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_slice.csv";
        final String centraPath = Resources.getResource("centra").getPath() + File.separator + "centra.csv";
        final String signaturePath = Resources.getResource("signature").getPath() + File.separator + "signature.png";
        final HmfReporterData reporterData =
                HmfReporterDataLoader.buildFromFiles(slicerPath, drupFilterPath, cosmicPath, centraPath, signaturePath);
        final JRDataSource dataSource = GenePanelDataSource.fromHmfReporterData(reporterData);
        assertNotNull(dataSource);
    }
}