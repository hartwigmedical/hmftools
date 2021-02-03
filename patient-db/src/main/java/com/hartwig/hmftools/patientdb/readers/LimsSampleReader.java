package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsAnalysisType;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientdb.data.ImmutableSampleData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LimsSampleReader {

    private static final Logger LOGGER = LogManager.getLogger(LimsSampleReader.class);

    @NotNull
    private final Lims lims;
    @NotNull
    private final Map<String, String> sampleToSetNameMap;
    @NotNull
    private final Set<String> sequencedSampleIds;

    public LimsSampleReader(@NotNull final Lims lims, @NotNull final Map<String, String> sampleToSetNameMap,
            @NotNull final Set<String> sequencedSampleIds) {
        this.lims = lims;
        this.sampleToSetNameMap = sampleToSetNameMap;
        this.sequencedSampleIds = sequencedSampleIds;
    }

    @Nullable
    public SampleData read(@NotNull String sampleBarcode, @NotNull String sampleId) {
        LocalDate arrivalDate = lims.arrivalDate(sampleBarcode, sampleId);
        boolean isSequenced = sequencedSampleIds.contains(sampleId);

        if (arrivalDate == null) {
            if (isSequenced) {
                LOGGER.warn("Could not find arrival date for sequenced sample {}", sampleId);
            }
            return null;
        }

        LocalDate samplingDate = lims.samplingDate(sampleBarcode);
        if (samplingDate == null && isSequenced && !lims.confirmedToHaveNoSamplingDate(sampleId)) {
            LOGGER.warn("Could not find sampling date for sequenced sample {}", sampleId);
        }

        String setName = sampleToSetNameMap.get(sampleId);
        if (setName == null && isSequenced) {
            LOGGER.warn("Could not resolve set name for sequenced sample {}", sampleId);
        }

        LimsCohortConfig cohortConfig = lims.cohortConfig(sampleBarcode);
        return ImmutableSampleData.builder()
                .sampleId(sampleId)
                .sampleBarcode(sampleBarcode)
                .cohortId(cohortConfig != null ? cohortConfig.cohortId() : Strings.EMPTY)
                .sequenced(isSequenced)
                .isSomaticTumorSample(lims.extractAnalysisType(sampleBarcode) == LimsAnalysisType.SOMATIC_T)
                .requiresCuratedPrimaryTumor(!lims.needsNoCuratedPrimaryTumor(lims.patientId(sampleBarcode)))
                .setName(setName != null ? setName : Strings.EMPTY)
                .arrivalDate(arrivalDate)
                .samplingDate(samplingDate)
                .dnaNanograms(lims.dnaNanograms(sampleBarcode))
                .limsPrimaryTumor(lims.primaryTumor(sampleBarcode))
                .pathologyTumorPercentage(lims.pathologyTumorPercentage(sampleBarcode))
                .pathologySampleId(lims.hospitalPathologySampleId(sampleBarcode))
                .build();
    }

    @Nullable
    public SampleData readSequencedSampleWithoutBarcode(@NotNull String sampleId) {
        LocalDate arrivalDate = lims.arrivalDate(Strings.EMPTY, sampleId);
        assert sequencedSampleIds.contains(sampleId);

        if (arrivalDate == null) {
            LOGGER.warn("Could not find arrival date for sequenced sample {}", sampleId);
            return null;
        }

        if (!lims.confirmedToHaveNoSamplingDate(sampleId)) {
            LOGGER.warn("Could not find sampling date for sequenced sample {}", sampleId);
        }

        String setName = sampleToSetNameMap.get(sampleId);
        if (setName == null) {
            LOGGER.warn("Could not resolve set name for sequenced sample {}", sampleId);
        }

        return ImmutableSampleData.builder()
                .sampleId(sampleId)
                .sampleBarcode(Strings.EMPTY)
                .cohortId(Strings.EMPTY)
                .sequenced(true)
                .isSomaticTumorSample(true)
                .requiresCuratedPrimaryTumor(true)
                .setName(setName != null ? setName : Strings.EMPTY)
                .arrivalDate(arrivalDate)
                .samplingDate(null)
                .dnaNanograms(null)
                .limsPrimaryTumor(null)
                .pathologyTumorPercentage(Lims.NOT_AVAILABLE_STRING)
                .pathologySampleId(Lims.NOT_AVAILABLE_STRING)
                .build();
    }
}
