package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.io.IOException;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class EcrfFileReader
{
    private EcrfFileReader()
    {
    }

    public static List<PatientData> readPatientData(@NotNull String pathToCsv) throws IOException
    {
        List<PatientData> result = new ArrayList<>();
        List<CsvEntry> entries = HmfCsvReader.read(pathToCsv);

        for(CsvEntry entry : entries)
        {
            var subjectKey = entry.get("SubjectKey").orElseThrow();
            var informedConsentDate = entry.get("DIC").map(LocalDate::parse);
            var registrationDate = entry.get("registrationDate").map(LocalDate::parse);
            var yearOfBirth = entry.get("YOB").map(Integer::parseInt);
            var gender = entry.get("GEN").map(i -> i.equals("0") ? PatientData.Gender.MALE : PatientData.Gender.FEMALE);
            var specimens = entry.get("CNS01").map(i -> i.equals("1"));
            var isDataAvailable = entry.get("CNS02").map(i -> i.equals("1"));
            var otherTrial = entry.get("OTHTRLC").map(i -> i.equals("Y"));
            var otherTrialCode = entry.get("DFHOTCC");
            var otherTrialDate = entry.get("DFOTSDTC").map(LocalDate::parse);
            var hadPreviousChemoTherapy = entry.get("CHYN").map(i -> i.equals("1"));
            var chemoEstimatedDate = entry.get("CHEDT");
            var chemoLastDoseDate = entry.get("CHDT").map(LocalDate::parse);
            var hadPreviousRadioTherapy = entry.get("RDYN").map(i -> i.equals("1"));
            var radioEstimatedDate = entry.get("RDEDT");
            var radioLastDoseDate = entry.get("RDDT").map(LocalDate::parse);

            var patientDataEntry = PatientData.builder()
                    .subjectKey(subjectKey)
                    .informedConsentDate(informedConsentDate)
                    .registrationDate(registrationDate)
                    .yearOfBirth(yearOfBirth)
                    .gender(gender)
                    .specimensCanBeStoredIndefinitelyAndBeUsedForFutureResearch(specimens)
                    .isDataAvailable(isDataAvailable)
                    .otherTrial(otherTrial)
                    .otherTrialCode(otherTrialCode)
                    .otherTrialDate(otherTrialDate)
                    .hadPreviousChemoTherapy(hadPreviousChemoTherapy)
                    .chemoEstimatedDate(chemoEstimatedDate)
                    .chemoLastDoseDate(chemoLastDoseDate)
                    .hadPreviousRadioTherapy(hadPreviousRadioTherapy)
                    .radioEstimatedDate(radioEstimatedDate)
                    .radioLastDoseDate(radioLastDoseDate)
                    .build();
            result.add(patientDataEntry);
        }
        return result;
    }

    private interface MappedColumn
    {
        String getMapping();
    }

    private enum BiospyDataColumn implements MappedColumn
    {
        COMBINED_KEY("Subject_FormRepeat_Keys"),
        SUBJECT_KEY("SubjectKey"),
        REGISTRATION_DATE("FormRepeatKey"),
        SAMPLE_DATE("BIOPTDTC"),
        SAMPLE_SITE("BIOPSITC"),
        SAMPLE_SITE_DETAILS("BIOPTISC"),
        SAMPLE_COLLECT_METHOD("SAMPTYPC"),
        STUDY_CODE("STCODEC"),
        OTHER_TRIAL("OTHTRLC"),
        OTHER_TRIAL_CODE("DFHOTCC"),
        OTHER_TRIAL_DATE("DFOTSDTC"),
        DIAGNOSIS("DFHDIAGC"),
        BDMWDPNR("BDMWDPNR"),
        T_NUMBER("BMDTNR"),
        WAS_WGS_SUCCESFUL("BMDBWSYN"),
        REASON_WGS_WAS_NOT_SUCCESFUL("BMDBWSNR"),
        SAMPLE_TYPE("BMDPTM"),
        WGS_REPORT_PIPELINE_VERSION("BMDWGSV"),
        HMF_REPORT_DATE("BMDWRPDT"),
        WGS_FINDINGS_SUMMARY("BMDWGSF"),
        TUMOR_TYPE_ONCO_TREE("BMDTTOT"),
        TUMOR_TYPE_HMF("BMDTTHMF");

        private final String csvColumnName;

        BiospyDataColumn(String csvColumnName)
        {
            this.csvColumnName = csvColumnName;
        }

        @Override
        public String getMapping()
        {
            return csvColumnName;
        }
    }

    private enum PrevTreatChemoColumn implements MappedColumn
    {
        SUBJECT_KEY("SubjectKey"),
        ITEM_GROUP_OID("ItemGroupOID"),
        ITEM_GROUP_REPEAT_KEY("ItemGroupRepeatKey"),
        MEDICAL_HISTORY_CATEGORY("CHCT"),
        CHEMO_DRUG_ONE("CHTR"),
        CHTS("CHTS"),
        CHEMO_DRUG_TWO("CHTR1"),
        CHTS1("CHTS1"),
        CHEMO_DRUG_THREE("CHTR2"),
        CHTS2("CHTS2"),
        CHEMO_DRUG_FOUR("CHTR3"),
        CHTS3("CHTS3");

        private final String csvColumnName;

        PrevTreatChemoColumn(String csvColumnName)
        {
            this.csvColumnName = csvColumnName;
        }

        @Override
        public String getMapping()
        {
            return csvColumnName;
        }
    }

    private enum PrevTreatRadColumn implements MappedColumn
    {
        SUBJECT_KEY("SubjectKey"),
        ITEM_GROUP_OID("ItemGroupOID"),
        ITEM_GROUP_REPEAT_KEY("ItemGroupRepeatKey"),
        RADIO_SITE("RDTR"),
        MEDICAL_HISTORY_CATEGORY("RDCT"),
        CUMULATIVE_DOSE("RDCD");

        private final String csvColumnName;

        PrevTreatRadColumn(String csvColumnName)
        {
            this.csvColumnName = csvColumnName;
        }

        @Override
        public String getMapping()
        {
            return csvColumnName;
        }
    }

    private enum TumorMeasureColumn implements MappedColumn
    {
        COMBINED_KEY("Subject_FormRepeat_Keys"),
        SUBJECT_KEY("SubjectKey"),
        FORM_REPEAT_KEY("FormRepeatKey"),
        MEASSURE_DATE("TMDTC"),
        TMDUMTXT("TMDUMTXT"),
        RECIST("RSPOA"),
        CONTINUE_TREATMENT("TMRESNR"),
        REASON_END_OF_TREATMENT("TMREOT"),
        REASON_END_OF_TREATMENT_SPECIFICATION("TMREOTSP");

        private final String csvColumnName;

        TumorMeasureColumn(String csvColumnName)
        {
            this.csvColumnName = csvColumnName;
        }

        @Override
        public String getMapping()
        {
            return csvColumnName;
        }
    }

    private enum TreatChemoAvlColumn implements MappedColumn
    {
        SUBJECT_KEY("SubjectKey"),
        CHEMO_CODE("TRAATCC"),
        TRASDT("TRASDT"),
        TRAEDT("TRAEDT");

        private final String csvColumnName;

        TreatChemoAvlColumn(String csvColumnName)
        {
            this.csvColumnName = csvColumnName;
        }

        @Override
        public String getMapping()
        {
            return csvColumnName;
        }
    }

}
