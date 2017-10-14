package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.ImmutableXMLEcrfDatamodel;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.lims.LimsJsonModel;

import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestUtil {

    public static final String SIGNATURE_PATH = Resources.getResource("signature").getPath() + File.separator + "signature.png";

    private PatientReporterTestUtil() {

    }

    @NotNull
    public static HmfReporterData testHmfReporterData() throws IOException, HartwigException {
        final String drupFilterPath = Resources.getResource("csv").getPath() + File.separator + "drup_genes.csv";
        final String cosmicPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_slice.csv";
        final String fusionPath = Resources.getResource("csv").getPath() + File.separator + "cosmic_gene_fusions.csv";
        return HmfReporterDataLoader.buildFromFiles(drupFilterPath, cosmicPath, fusionPath);
    }

    @NotNull
    public static BaseReporterData testBaseReporterData() throws IOException, EmptyFileException {
        final String centerPath = Resources.getResource("center").getPath() + File.separator + "centers.csv";
        final CpctEcrfModel ecrfModel = new CpctEcrfModel(
                ImmutableXMLEcrfDatamodel.of(Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList(),
                        Lists.newArrayList()), Lists.newArrayList());
        final LimsJsonModel limsModel = LimsJsonModel.buildEmptyModel();
        final CenterModel centerModel = Center.readFromCSV(centerPath);

        return ImmutableBaseReporterData.of(ecrfModel, limsModel, centerModel, SIGNATURE_PATH);
    }
}
