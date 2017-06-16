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
        final String slicerPath = Resources.getResource("bed").getPath() + File.separator + "HMF_Slicing_v2.tsv";
        final String drupFilterPath = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
        final String cosmicPath = Resources.getResource("csv").getPath() + File.separator + "cosmic.csv";
        final HmfReporterData reporterData = HmfReporterDataLoader.buildFromFiles(slicerPath, drupFilterPath,
                cosmicPath);
        final JRDataSource dataSource = GenePanelDataSource.fromHmfReporterData(reporterData);
        assertNotNull(dataSource);
    }
}