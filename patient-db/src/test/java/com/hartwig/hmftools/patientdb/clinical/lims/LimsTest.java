package com.hartwig.hmftools.patientdb.clinical.lims;

import static com.hartwig.hmftools.patientdb.clinical.lims.LimsTestUtil.createLimsSampleDataBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.patientdb.clinical.lims.cohort.LimsCohortTestFactory;
import com.hartwig.hmftools.patientdb.clinical.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.patientdb.clinical.lims.hospital.HospitalModel;
import com.hartwig.hmftools.patientdb.clinical.lims.hospital.ImmutableHospitalAddress;
import com.hartwig.hmftools.patientdb.clinical.lims.hospital.ImmutableHospitalModel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LimsTest {

    private static final String REF_SAMPLE_BARCODE = "FR1234";
    private static final String TUMOR_SAMPLE_BARCODE = "FR1";
    private static final String TUMOR_SAMPLE_ID = "CPCT-HMF-Test-Sample";
    private static final String TEST_COHORT = "CPCT";
    private static final String SUBMISSION = "ABCDEF123";

    @Test
    public void canReadProperlyDefinedSample() {
        String patientId = "CPCT02991111";
        String arrivalDate = "2017-05-01";
        String samplingDate = "2017-04-15";
        String dnaConcentration = "10";
        String purityShallowSeq = "0.71";
        String primaryTumor = "Prostate";
        String biopsyLocation = "Skin";
        String labSopVersions = "PREP1V2-QC1V2-SEQ1V2";
        String cohort = "CPCT";
        String projectName = "projectX";
        String requesterName = "henk";
        String requesterEmail = "henk@hmf.nl";
        String hospitalPatientId = "Henkie";
        String hospitalPathologySampleId = "Henkie's sample";
        boolean reportGermlineVariants = false;
        boolean reportViralPresence = false;

        LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(TUMOR_SAMPLE_ID)
                .patientId(patientId)
                .arrivalDate(arrivalDate)
                .samplingDate(samplingDate)
                .dnaConcentration(dnaConcentration)
                .primaryTumor(primaryTumor)
                .biopsySite(biopsyLocation)
                .labSopVersions(labSopVersions)
                .submission(SUBMISSION)
                .cohort(cohort)
                .refBarcode(REF_SAMPLE_BARCODE)
                .tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .shallowSeq(true)
                .hospitalPatientId(hospitalPatientId)
                .hospitalPathologySampleId(hospitalPathologySampleId)
                .germlineReportingLevel(Strings.EMPTY)
                .reportGermlineVariants(reportGermlineVariants)
                .reportViralPresence(reportViralPresence)
                .build();

        LimsJsonSubmissionData submissionData = ImmutableLimsJsonSubmissionData.builder()
                .submission(SUBMISSION)
                .projectName(projectName)
                .reportContactName(requesterName)
                .reportContactEmail(requesterEmail)
                .build();

        LimsShallowSeqData shallowSeqData = ImmutableLimsShallowSeqData.builder()
                .sampleBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId(TUMOR_SAMPLE_ID)
                .purityShallowSeq(purityShallowSeq)
                .hasReliableQuality(true)
                .hasReliablePurity(true)
                .build();

        Lims lims = buildFullTestLims(sampleData, submissionData, shallowSeqData);

        assertEquals(1, lims.sampleBarcodeCount());
        assertEquals(1, lims.sampleBarcodes().size());
        assertFalse(lims.confirmedToHaveNoSamplingDate(TUMOR_SAMPLE_ID));

        lims.validateSampleBarcodeCombination(REF_SAMPLE_BARCODE, Strings.EMPTY, TUMOR_SAMPLE_BARCODE, TUMOR_SAMPLE_ID);

        assertEquals(patientId, lims.patientId(TUMOR_SAMPLE_BARCODE));
        assertEquals(TUMOR_SAMPLE_ID, lims.sampleId(TUMOR_SAMPLE_BARCODE));
        assertEquals(LimsTestUtil.toDate(arrivalDate), lims.arrivalDate(TUMOR_SAMPLE_BARCODE, TUMOR_SAMPLE_ID));
        assertEquals(LimsTestUtil.toDate(samplingDate), lims.samplingDate(TUMOR_SAMPLE_BARCODE));
        assertFalse(lims.isBlacklisted(patientId));

        assertEquals(SUBMISSION, lims.submissionId(TUMOR_SAMPLE_BARCODE));
        assertEquals(projectName, lims.projectName(TUMOR_SAMPLE_BARCODE));
        assertEquals(requesterName, lims.requesterName(TUMOR_SAMPLE_BARCODE));
        assertEquals(requesterEmail, lims.requesterEmail(TUMOR_SAMPLE_BARCODE));

        Integer dnaAmount = lims.dnaNanograms(TUMOR_SAMPLE_BARCODE);
        assertNotNull(dnaAmount);
        assertEquals(500L, (int) dnaAmount);

        assertEquals("71%", lims.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.pathologyTumorPercentage(TUMOR_SAMPLE_BARCODE));
        assertEquals(primaryTumor, lims.primaryTumor(TUMOR_SAMPLE_BARCODE));
        assertEquals(biopsyLocation, lims.biopsyLocation(TUMOR_SAMPLE_BARCODE));
        assertEquals(labSopVersions, lims.labProcedures(TUMOR_SAMPLE_BARCODE));

        assertEquals(hospitalPatientId, lims.hospitalPatientId(TUMOR_SAMPLE_BARCODE));
        assertEquals(hospitalPathologySampleId, lims.hospitalPathologySampleId(TUMOR_SAMPLE_BARCODE));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING, lims.germlineReportingChoice(TUMOR_SAMPLE_BARCODE, false));

        assertEquals(LimsGermlineReportingLevel.NO_REPORTING, lims.germlineReportingChoice(TUMOR_SAMPLE_BARCODE, false));
        assertEquals(reportGermlineVariants, lims.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
        assertEquals(reportViralPresence, lims.reportViralPresence(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void worksForNonExistingSamplesAndSubmissions() {
        Lims lims = LimsFactory.empty();
        String doesNotExistSample = "DoesNotExist";

        lims.validateSampleBarcodeCombination(doesNotExistSample, doesNotExistSample, doesNotExistSample, doesNotExistSample);
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.patientId(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.sampleId(doesNotExistSample));
        assertNull(lims.arrivalDate(doesNotExistSample, doesNotExistSample));
        assertNull(lims.samplingDate(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.submissionId(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.projectName(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.requesterName(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.requesterEmail(doesNotExistSample));
        assertNull(lims.dnaNanograms(doesNotExistSample));
        assertEquals(Lims.NOT_PERFORMED_STRING, lims.purityShallowSeq(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.pathologyTumorPercentage(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.primaryTumor(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.labProcedures(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.hospitalPatientId(doesNotExistSample));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.hospitalPathologySampleId(doesNotExistSample));
        assertEquals(LimsGermlineReportingLevel.NO_REPORTING, lims.germlineReportingChoice(doesNotExistSample, false));
        assertFalse(lims.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
        assertFalse(lims.reportViralPresence(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void fallBackOnPreLimsArrivalDateWorks() {
        LocalDate date = LimsTestUtil.toDate("2017-10-03");

        Lims lims = buildTestLimsWithPreLimsArrivalDateForSampleId(TUMOR_SAMPLE_ID, date);

        assertEquals(date, lims.arrivalDate("DoesNotExist", TUMOR_SAMPLE_ID));
    }

    @Test
    public void invalidDataYieldsNullOrNA() {
        LimsJsonSampleData sampleData = createLimsSampleDataBuilder().sampleId(TUMOR_SAMPLE_ID)
                .tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .arrivalDate("IsNotADate")
                .samplingDate(null)
                .dnaConcentration("IsNotADNAConcentration")
                .pathologyTumorPercentage("IsNotANumber")
                .labSopVersions("anything")
                .reportViralPresence(false)
                .cohort("CPCT")
                .build();

        Lims lims = buildTestLimsWithSample(sampleData);

        assertEquals(1, lims.sampleBarcodeCount());

        assertNull(lims.arrivalDate(TUMOR_SAMPLE_BARCODE, TUMOR_SAMPLE_ID));
        assertNull(lims.samplingDate(TUMOR_SAMPLE_BARCODE));
        assertNull(lims.dnaNanograms(TUMOR_SAMPLE_BARCODE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.pathologyTumorPercentage(TUMOR_SAMPLE_BARCODE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims.labProcedures(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canReadHospitalData() {
        LimsJsonSampleData sample =
                createLimsSampleDataBuilder().sampleId(TUMOR_SAMPLE_ID).tumorBarcode(TUMOR_SAMPLE_BARCODE).cohort(TEST_COHORT).build();

        Lims emptyHospitalModel = buildTestLimsWithWithHospitalModel(sample, ImmutableHospitalModel.builder().build());

        HospitalContactData emptyContact = emptyHospitalModel.hospitalContactData(TUMOR_SAMPLE_BARCODE);
        assertEquals(Lims.NOT_AVAILABLE_STRING, emptyContact.hospitalPI());
        assertEquals(Lims.NOT_AVAILABLE_STRING, emptyContact.requesterName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, emptyContact.requesterEmail());
        assertEquals(Lims.NOT_AVAILABLE_STRING, emptyContact.hospitalName());
        assertEquals(Lims.NOT_AVAILABLE_STRING, emptyContact.hospitalAddress());

        HospitalModel testHospitalModel = ImmutableHospitalModel.builder()
                .putHospitalAddressMap("1", ImmutableHospitalAddress.of("Name", "Zip", "City"))
                .putSampleToHospitalMapping(TUMOR_SAMPLE_ID, "1")
                .build();

        Lims withHospitalModel = buildTestLimsWithWithHospitalModel(sample, testHospitalModel);
        HospitalContactData contact = withHospitalModel.hospitalContactData(TUMOR_SAMPLE_BARCODE);
        assertEquals("Name", contact.hospitalName());
        assertEquals("Zip City", contact.hospitalAddress());
    }

    @Test
    public void canExtractLimsViralInsertionsChoice() {
        LimsJsonSampleData sampleDataTrue = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId(TUMOR_SAMPLE_ID)
                .cohort(TEST_COHORT)
                .reportViralPresence(true)
                .build();
        LimsJsonSampleData sampleDataFalse = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId(TUMOR_SAMPLE_ID)
                .cohort(TEST_COHORT)
                .reportViralPresence(false)
                .build();

        Lims limsTrue = buildTestLimsWithSample(sampleDataTrue);
        Lims limsFalse = buildTestLimsWithSample(sampleDataFalse);

        assertTrue(limsTrue.reportViralPresence(TUMOR_SAMPLE_BARCODE));

        assertFalse(limsFalse.reportViralPresence(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canExtractLimsViralInsertionsChoiceException() {
        LimsJsonSampleData sampleDataTrue = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId(TUMOR_SAMPLE_ID)
                .cohort(TEST_COHORT)
                .reportViralPresence(true)
                .build();
        Lims limsTrue = buildTestLimsWithSample(sampleDataTrue);

        assertFalse(limsTrue.reportViralPresence("does not exist"));
    }

    @Test
    public void canExtractLimsReportableGermlineVariants() {
        LimsJsonSampleData sampleDataCPCTTrue = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId("CPCT02991111T")
                .reportGermlineVariants(true)
                .cohort("CPCT")
                .build();
        LimsJsonSampleData sampleDataCPCTFalse = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId("CPCT02991111T")
                .reportGermlineVariants(false)
                .cohort("CPCT")
                .build();

        LimsJsonSampleData sampleDataWIDETrue = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId("WIDE01010000T")
                .reportGermlineVariants(true)
                .cohort("WIDE")
                .build();
        LimsJsonSampleData sampleDataWIDEFalse = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId("WIDE01010000T")
                .reportGermlineVariants(false)
                .cohort("WIDE")
                .build();

        Lims limsCPCTTrue = buildTestLimsWithSample(sampleDataCPCTTrue);
        Lims limsCPCTFalse = buildTestLimsWithSample(sampleDataCPCTFalse);
        Lims limsWIDETrue = buildTestLimsWithSample(sampleDataWIDETrue);
        Lims limsWIDEFalse = buildTestLimsWithSample(sampleDataWIDEFalse);

        assertTrue(limsCPCTTrue.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
        assertFalse(limsCPCTFalse.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
        assertTrue(limsWIDETrue.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
        assertFalse(limsWIDEFalse.reportGermlineVariants(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canExtractLimsReportableGermlineVariantsException() {
        LimsJsonSampleData sampleDataCPCTTrue = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .sampleId("CPCT02991111T")
                .reportGermlineVariants(true)
                .cohort("CPCT")
                .build();
        Lims limsCPCTTrue = buildTestLimsWithSample(sampleDataCPCTTrue);

        assertFalse(limsCPCTTrue.reportGermlineVariants("does not exist"));
    }

    @Test
    public void missingOrMalformedShallowSeqDataForSampleYieldsNA() {
        LimsJsonSampleData sampleData1 = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE).shallowSeq(true).build();

        Lims lims1 = buildTestLimsWithSample(sampleData1);
        assertEquals(Lims.NOT_PERFORMED_STRING, lims1.purityShallowSeq(TUMOR_SAMPLE_BARCODE));

        LimsJsonSampleData sampleData2 = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE).shallowSeq(true).build();
        Lims lims2 = buildTestLimsWithSampleAndShallowSeq(sampleData2, "NotANumber", true, true);
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims2.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canRetrievePathologyPercentageForSample() {
        LimsJsonSampleData sampleData1 =
                createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE).shallowSeq(false).pathologyTumorPercentage("70").build();

        Lims lims1 = buildTestLimsWithSample(sampleData1);
        assertEquals(Lims.NOT_PERFORMED_STRING, lims1.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
        assertEquals("70%", lims1.pathologyTumorPercentage(TUMOR_SAMPLE_BARCODE));

        LimsJsonSampleData sampleData2 = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE)
                .shallowSeq(false)
                .pathologyTumorPercentage("NotANumber")
                .build();

        Lims lims2 = buildTestLimsWithSample(sampleData2);
        assertEquals(Lims.NOT_PERFORMED_STRING, lims2.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
        assertEquals(Lims.NOT_AVAILABLE_STRING, lims2.pathologyTumorPercentage(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canRetrieveShallowSeqPurityForSample() {
        LimsJsonSampleData sampleData = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE).shallowSeq(true).build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "0.2", false, true);

        assertEquals("20%", lims.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
    }

    @Test
    public void canRetrieveShallowSeqBelowDetectionLimitForSample() {
        LimsJsonSampleData sampleData = createLimsSampleDataBuilder().tumorBarcode(TUMOR_SAMPLE_BARCODE).shallowSeq(true).build();

        Lims lims = buildTestLimsWithSampleAndShallowSeq(sampleData, "0.10", true, false);
        assertEquals(Lims.PURITY_NOT_RELIABLE_STRING, lims.purityShallowSeq(TUMOR_SAMPLE_BARCODE));
    }

    @NotNull
    private static Lims buildFullTestLims(@NotNull LimsJsonSampleData sample, @NotNull LimsJsonSubmissionData submissionData,
            @NotNull LimsShallowSeqData shallowSeqData) {
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = Maps.newHashMap();
        dataPerSampleBarcode.put(sample.tumorBarcode(), sample);

        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        dataPerSubmission.put(submissionData.submission(), submissionData);

        Map<String, LimsShallowSeqData> shallowSeqDataPerSampleBarcode = Maps.newHashMap();
        shallowSeqDataPerSampleBarcode.put(shallowSeqData.sampleBarcode(), shallowSeqData);

        Map<String, LocalDate> preLimsArrivalDatesPerSampleId = Maps.newHashMap();
        Set<String> sampleIdsWithoutSamplingDates = Sets.newHashSet();
        Set<String> blacklistedPatients = Sets.newHashSet();
        Set<String> patientsWithoutCuratedPrimaryTumor = Sets.newHashSet();
        HospitalModel hospitalModel = ImmutableHospitalModel.builder().build();

        return new Lims(buildCohortModelFromSampleCohort(sample.cohort()),
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqDataPerSampleBarcode,
                preLimsArrivalDatesPerSampleId,
                sampleIdsWithoutSamplingDates,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    private static Lims buildTestLimsWithSample(@NotNull LimsJsonSampleData sample) {
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = Maps.newHashMap();
        dataPerSampleBarcode.put(sample.tumorBarcode(), sample);

        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSampleBarcode = Maps.newHashMap();
        Map<String, LocalDate> preLimsArrivalDatesPerSampleId = Maps.newHashMap();
        Set<String> sampleIdsWithoutSamplingDate = Sets.newHashSet();
        Set<String> blacklistedPatients = Sets.newHashSet();
        Set<String> patientsWithoutCuratedPrimaryTumor = Sets.newHashSet();
        HospitalModel hospitalModel = ImmutableHospitalModel.builder().build();

        return new Lims(buildCohortModelFromSampleCohort(sample.cohort()),
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqDataPerSampleBarcode,
                preLimsArrivalDatesPerSampleId,
                sampleIdsWithoutSamplingDate,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    private static Lims buildTestLimsWithSampleAndShallowSeq(@NotNull LimsJsonSampleData sample, @NotNull String shallowSeqPurity,
            boolean hasReliableQuality, boolean hasReliablePurity) {
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = Maps.newHashMap();
        dataPerSampleBarcode.put(sample.tumorBarcode(), sample);

        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSampleBarcode = Maps.newHashMap();
        shallowSeqDataPerSampleBarcode.put(sample.tumorBarcode(),
                ImmutableLimsShallowSeqData.of("FR1", sample.sampleId(), shallowSeqPurity, hasReliableQuality, hasReliablePurity));

        Map<String, LocalDate> preLimsArrivalDatesPerSampleId = Maps.newHashMap();
        Set<String> sampleIdsWithoutSamplingDate = Sets.newHashSet();
        Set<String> blacklistedPatients = Sets.newHashSet();
        Set<String> patientsWithoutCuratedPrimaryTumor = Sets.newHashSet();
        HospitalModel hospitalModel = ImmutableHospitalModel.builder().build();

        return new Lims(buildCohortModelFromSampleCohort(sample.cohort()),
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqDataPerSampleBarcode,
                preLimsArrivalDatesPerSampleId,
                sampleIdsWithoutSamplingDate,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    private static Lims buildTestLimsWithPreLimsArrivalDateForSampleId(@NotNull String sampleId, @NotNull LocalDate date) {
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = Maps.newHashMap();
        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSampleBarcode = Maps.newHashMap();

        Map<String, LocalDate> preLimsArrivalDatesPerSampleId = Maps.newHashMap();
        preLimsArrivalDatesPerSampleId.put(sampleId, date);

        Set<String> sampleIdsWithoutSamplingDate = Sets.newHashSet();
        Set<String> blacklistedPatients = Sets.newHashSet();
        Set<String> patientsWithoutCuratedPrimaryTumor = Sets.newHashSet();
        HospitalModel hospitalModel = ImmutableHospitalModel.builder().build();
        LimsCohortModel cohortModel = ImmutableLimsCohortModel.builder().build();

        return new Lims(cohortModel,
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqDataPerSampleBarcode,
                preLimsArrivalDatesPerSampleId,
                sampleIdsWithoutSamplingDate,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    private static Lims buildTestLimsWithWithHospitalModel(@NotNull LimsJsonSampleData sample, @NotNull HospitalModel hospitalModel) {
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = Maps.newHashMap();
        dataPerSampleBarcode.put(sample.tumorBarcode(), sample);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = Maps.newHashMap();
        Map<String, LimsShallowSeqData> shallowSeqDataPerSampleBarcode = Maps.newHashMap();

        Map<String, LocalDate> preLimsArrivalDatesPerSampleId = Maps.newHashMap();
        Set<String> sampleIdsWithoutSamplingDate = Sets.newHashSet();
        Set<String> blacklistedPatients = Sets.newHashSet();
        Set<String> patientsWithoutCuratedPrimaryTumor = Sets.newHashSet();

        return new Lims(buildCohortModelFromSampleCohort(sample.cohort()),
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqDataPerSampleBarcode,
                preLimsArrivalDatesPerSampleId,
                sampleIdsWithoutSamplingDate,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    private static LimsCohortModel buildCohortModelFromSampleCohort(@NotNull String sampleCohort) {
        Map<String, LimsCohortConfig> configMap = Maps.newHashMap();
        configMap.put(sampleCohort, LimsCohortTestFactory.createAllDisabledCohortConfig(sampleCohort));
        return ImmutableLimsCohortModel.builder().limsCohortMap(configMap).build();
    }
}