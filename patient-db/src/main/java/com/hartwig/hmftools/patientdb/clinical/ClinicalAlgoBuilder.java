package com.hartwig.hmftools.patientdb.clinical;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatusReader;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfFileReader;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideFiveDays;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WidePreAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideResponseData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ClinicalAlgoBuilder {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalAlgoBuilder.class);

    private ClinicalAlgoBuilder() {
    }

    @NotNull
    public static ClinicalAlgo fromConfig(@NotNull ClinicalAlgoConfig config) throws IOException, XMLStreamException {
        List<DoidNode> doidNodes = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJson()).nodes();
        PrimaryTumorCurator primaryTumorCurator =
                new PrimaryTumorCurator(config.tumorLocationMappingTsv(), config.tumorLocationOverridesTsv(), doidNodes);
        BiopsySiteCurator biopsySiteCurator = new BiopsySiteCurator(config.biopsyMappingTsv());
        TreatmentCurator treatmentCurator = new TreatmentCurator(config.treatmentMappingTsv());

        EcrfModels ecrfModels = loadEcrfModels(config);

        return new ClinicalAlgo(ecrfModels, primaryTumorCurator, biopsySiteCurator, treatmentCurator);
    }

    @NotNull
    private static EcrfModels loadEcrfModels(@NotNull ClinicalAlgoConfig config) throws IOException, XMLStreamException {
        EcrfModel cpctEcrfModel = buildCpctEcrfModel(config);
        EcrfModel drupEcrfModel = buildDrupEcrfModel(config);
        WideEcrfModel wideEcrfModel = buildWideEcrfModel(config);

        return ImmutableEcrfModels.builder().cpctModel(cpctEcrfModel).drupModel(drupEcrfModel).wideModel(wideEcrfModel).build();
    }

    @NotNull
    private static EcrfModel buildCpctEcrfModel(@NotNull ClinicalAlgoConfig config) throws IOException, XMLStreamException {
        LOGGER.info("Loading CPCT eCRF from {}", config.cpctEcrfFile());
        FormStatusModel cpctFormStatusModel = FormStatusReader.buildModelFromCsv(config.cpctFormStatusCsv());
        EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLWithFormStates(config.cpctEcrfFile(), cpctFormStatusModel);
        LOGGER.info(" Finished loading CPCT eCRF. Read {} patients", cpctEcrfModel.patientCount());

        return cpctEcrfModel;
    }

    @NotNull
    private static EcrfModel buildDrupEcrfModel(@NotNull ClinicalAlgoConfig config) throws FileNotFoundException, XMLStreamException {
        LOGGER.info("Loading DRUP eCRF from {}", config.drupEcrfFile());
        EcrfModel drupEcrfModel = EcrfModel.loadFromXMLNoFormStates(config.drupEcrfFile());
        LOGGER.info(" Finished loading DRUP eCRF. Read {} patients", drupEcrfModel.patientCount());

        return drupEcrfModel;
    }

    @NotNull
    private static WideEcrfModel buildWideEcrfModel(@NotNull ClinicalAlgoConfig config) throws IOException {
        WideEcrfModel wideEcrfModel;

        if (config.doProcessWideClinicalData()) {
            LOGGER.info("Loading WIDE eCRF");

            String preAvlTreatmentCsv = config.widePreAvlTreatmentCsv();
            List<WidePreAvlTreatmentData> preAvlTreatments = WideEcrfFileReader.readPreAvlTreatments(preAvlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE pre-AVL-treatments from {}", preAvlTreatments.size(), preAvlTreatmentCsv);

            String biopsyCsv = config.wideBiopsyCsv();
            List<WideBiopsyData> biopsies = WideEcrfFileReader.readBiopsies(biopsyCsv);
            LOGGER.info(" Loaded {} WIDE biopsies from {}", biopsies.size(), biopsyCsv);

            String avlTreatmentCsv = config.wideAvlTreatmentCsv();
            List<WideAvlTreatmentData> avlTreatments = WideEcrfFileReader.readAvlTreatments(avlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE AVL treatments from {}", avlTreatments.size(), avlTreatmentCsv);

            String wideResponseCsv = config.wideResponseCsv();
            List<WideResponseData> responses = WideEcrfFileReader.readResponses(wideResponseCsv);
            LOGGER.info(" Loaded {} WIDE responses from {}", responses.size(), wideResponseCsv);

            String fiveDaysCsv = config.wideFiveDaysCsv();
            List<WideFiveDays> fiveDays = WideEcrfFileReader.readFiveDays(fiveDaysCsv);
            LOGGER.info(" Loaded {} WIDE five days entries from {}", fiveDays.size(), fiveDaysCsv);

            wideEcrfModel = ImmutableWideEcrfModel.builder()
                    .preAvlTreatments(preAvlTreatments)
                    .biopsies(biopsies)
                    .avlTreatments(avlTreatments)
                    .responses(responses)
                    .fiveDays(fiveDays)
                    .build();
        } else {
            LOGGER.info("Skipping the loading of WIDE eCRF");
            wideEcrfModel = ImmutableWideEcrfModel.builder()
                    .preAvlTreatments(Lists.newArrayList())
                    .biopsies(Lists.newArrayList())
                    .avlTreatments(Lists.newArrayList())
                    .responses(Lists.newArrayList())
                    .fiveDays(Lists.newArrayList())
                    .build();
        }

        return wideEcrfModel;
    }
}
