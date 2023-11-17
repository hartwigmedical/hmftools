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
            // required fields
            var subjectKey = entry.get("SubjectKey").orElseThrow();

            // optional fields
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

    public static List<BiopsyData> readBiopsyData(@NotNull String pathToCsv) throws IOException
    {
        List<BiopsyData> result = new ArrayList<>();

        List<CsvEntry> entries = HmfCsvReader.read(pathToCsv);

        for(CsvEntry entry : entries)
        {
            // required fields
            var combinedKey = entry.get("Subject_FormRepeat_Keys").orElseThrow();
            var subjectKey = entry.get("SubjectKey").orElseThrow();
            var registrationDate = entry.get("FormRepeatKey").map(LocalDate::parse).orElseThrow();

            // optional fields
            var sampleDate = entry.get("BIOPTDTC").map(LocalDate::parse);
            var sampleSite = entry.get("BIOPSITC");
            var sampleSiteDetails = entry.get("BIOPTISC");
            var sampleCollectionMethod = entry.get("SAMPTYPC");
            var studyCode = entry.get("STCODEC");
            var otherTrial = entry.get("OTHTRLC");
            var otherTrialCode = entry.get("DFHOTCC");
            var otherTrialDate = entry.get("DFOTSDTC").map(LocalDate::parse);
            var diagnosis = entry.get("DFHDIAGC");
            var BDMWDPNR = entry.get("BDMWDPNR");
            var tNumber = entry.get("BMDTNR");
            var wasWgsSuccessful = entry.get("BMDBWSYN").map(i -> i.equals("1"));
            var reasonWgsWasNotSuccessful = entry.get("BMDBWSNR");
            var sampleType = entry.get("BMDPTM");
            var wgsReportPipelineVersion = entry.get("BMDWGSV");
            var hmfReportDate = entry.get("BMDWRPDT").map(LocalDate::parse);
            var wgsFindingsSummary = entry.get("BMDWGSF");
            var tumorTypeOncoTree = entry.get("BMDTTOT");
            var tumorTypeHmf = entry.get("BMDTTHMF");

            var biopsyDataEntry = BiopsyData.builder()
                    .combinedKey(combinedKey)
                    .subjectKey(subjectKey)
                    .registrationDate(registrationDate)
                    .sampleDate(sampleDate)
                    .sampleSite(sampleSite)
                    .sampleSiteDetails(sampleSiteDetails)
                    .sampleCollectMethod(sampleCollectionMethod)
                    .studyCode(studyCode)
                    .otherTrial(otherTrial)
                    .otherTrialCode(otherTrialCode)
                    .otherTrialDate(otherTrialDate)
                    .diagnosis(diagnosis)
                    .BDMWDPNR(BDMWDPNR)
                    .tNumber(tNumber)
                    .wasWgsSuccessful(wasWgsSuccessful)
                    .reasonWgsWasNotSuccessful(reasonWgsWasNotSuccessful)
                    .sampleType(sampleType)
                    .wgsReportPipelineVersion(wgsReportPipelineVersion)
                    .hmfReportDate(hmfReportDate)
                    .wgsFindingsSummary(wgsFindingsSummary)
                    .tumorTypeOncoTree(tumorTypeOncoTree)
                    .tumorTypeHmf(tumorTypeHmf)
                    .build();

            result.add(biopsyDataEntry);
        }
        return result;
    }

    public static List<PrevTreatChemoData> readPrevTreatChemoData(@NotNull String pathToCsv) throws IOException
    {
        List<PrevTreatChemoData> result = new ArrayList<>();
        List<CsvEntry> entries = HmfCsvReader.read(pathToCsv);

        for(CsvEntry entry : entries)
        {
            // required fields
            var subjectKey = entry.get("SubjectKey").orElseThrow();
            var itemGroupOid = entry.get("ItemGroupOID").orElseThrow();
            var itemGroupRepeatKey = entry.get("ItemGroupRepeatKey").map(Integer::parseInt).orElseThrow();

            // optional fields
            var medicalHistoryCategory = entry.get("CHCT");
            var chemoDrugOne = entry.get("CHTR").map(Integer::parseInt);
            var CHTS = entry.get("CHTS").map(Integer::parseInt);
            var chemoDrugTwo = entry.get("CHTR1").map(Integer::parseInt);
            var CHTS1 = entry.get("CHTS1").map(Integer::parseInt);
            var chemoDrugThree = entry.get("CHTR2").map(Integer::parseInt);
            var CHTS2 = entry.get("CHTS2").map(Integer::parseInt);
            var chemoDrugFour = entry.get("CHTR3").map(Integer::parseInt);
            var CHTS3 = entry.get("CHTS3").map(Integer::parseInt);

            var prevTreatChemoDataEntry = PrevTreatChemoData.builder()
                    .subjectKey(subjectKey)
                    .itemGroupOid(itemGroupOid)
                    .itemGroupRepeatedKey(itemGroupRepeatKey)
                    .medicalHistoryCategory(medicalHistoryCategory)
                    .chemoDrugOne(chemoDrugOne)
                    .CHTS(CHTS)
                    .chemoDrugTwo(chemoDrugTwo)
                    .CHTS1(CHTS1)
                    .chemoDrugThree(chemoDrugThree)
                    .CHTS2(CHTS2)
                    .chemoDrugFour(chemoDrugFour)
                    .CHTS3(CHTS3)
                    .build();
            result.add(prevTreatChemoDataEntry);
        }
        return result;
    }

    public static List<PrevTreatRadData> readPrevTreatRadData(@NotNull String pathToCsv) throws IOException
    {
        List<PrevTreatRadData> result = new ArrayList<>();

        List<CsvEntry> entries = HmfCsvReader.read(pathToCsv);

        for(CsvEntry entry : entries)
        {
            // required fields
            var subjectKey = entry.get("SubjectKey").orElseThrow();
            var itemGroupOid = entry.get("ItemGroupOID").orElseThrow();
            var itemGroupRepeatKey = entry.get("ItemGroupRepeatKey").map(Integer::parseInt).orElseThrow();

            // optional fields
            var radioSite = entry.get("RDTR");
            var medicalHistoryCategory = entry.get("RDCT");
            var cumulativeDose = entry.get("RDCD").map(Integer::parseInt);

            var prevTreatRadDataEntry = PrevTreatRadData.builder()
                    .subjectKey(subjectKey)
                    .itemGroupOid(itemGroupOid)
                    .itemGroupRepeatKey(itemGroupRepeatKey)
                    .radioSite(radioSite)
                    .medicalHistoryCategory(medicalHistoryCategory)
                    .cumulativeDose(cumulativeDose)
                    .build();

            result.add(prevTreatRadDataEntry);
        }
        return result;
    }

    public static List<TumorMeasureData> readTumorMeasureData(@NotNull String pathToCsv) throws IOException
    {
        List<TumorMeasureData> result = new ArrayList<>();
        List<CsvEntry> entries = HmfCsvReader.read(pathToCsv);

        for(CsvEntry entry : entries)
        {
            // required fields
            var combinedKey = entry.get("Subejct_FormRepeat_Keys").orElseThrow();
            var subjectKey = entry.get("SubjectKey").orElseThrow();
            var formRepeatKey = entry.get("FormRepeatKey").orElseThrow();

            // optional fields
            var measureDate = entry.get("TMDTC").map(LocalDate::parse);
            var TMDUMTXT = entry.get("TMDUMTXT");
            var recist = entry.get("RSPOA").map(Integer::parseInt);
            var continueTreatment = entry.get("TMREOT").map(i -> i.equals("1"));
            var reasonEndOfTreatment = entry.get("TMREOT").map(TumorMeasureData.EndOfTreatmentReason::valueOf);
            var reasonEndOfTreatmentSpecification = entry.get("TMREOTSP");

            var tumorMeasureDataEntry = TumorMeasureData.builder()
                    .combinedKey(combinedKey)
                    .subjectKey(subjectKey)
                    .formRepeatKey(formRepeatKey)
                    .measureDate(measureDate)
                    .TMDUMTXT(TMDUMTXT)
                    .recist(recist)
                    .continueTreatment(continueTreatment)
                    .reasonEndOfTreatment(reasonEndOfTreatment)
                    .reasonEndOfTreatmentSpecification(reasonEndOfTreatmentSpecification)
                    .build();

            result.add(tumorMeasureDataEntry);
        }
        return result;
    }

    private interface MappedColumn
    {
        String getMapping();
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
