package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.PredictionData.DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class DataLoader
{
    public static final String LILAC_COVERAGE_FILE_ID = ".lilac.csv";
    public static final String MHCF_PREDICTION_FILE_ID = ".mcf.predictions.csv";

    public static List<AlleleCoverage> loadAlleleCoverage(final String sampleId, final String lilacDir)
    {
        List<AlleleCoverage> alleleCoverages = Lists.newArrayListWithExpectedSize(EXPECTED_ALLELE_COUNT);

        try
        {
            String filename = lilacDir + sampleId + LILAC_COVERAGE_FILE_ID;

            List<String> lines = Files.readAllLines(new File(filename).toPath());

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), AlleleCoverage.DELIM);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get("Allele");
            int tumorCnIndex = fieldsIndexMap.get("TumorCopyNumber");

            // only select the somatic mutations which are predicted to silence/disable the allele
            List<Integer> somVariantIndices = Lists.newArrayList(
                    fieldsIndexMap.get("SomaticMissense"),
                    fieldsIndexMap.get("SomaticSplice"),
                    fieldsIndexMap.get("SomaticNonsenseOrFrameshift"),
                    fieldsIndexMap.get("SomaticInframeIndel"));

            for(String line : lines)
            {
                alleleCoverages.add(AlleleCoverage.fromCsv(line, alleleIndex, tumorCnIndex, somVariantIndices));
            }
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) Lilac allele coverage file: {}", sampleId, exception.toString());
        }

        return alleleCoverages;
    }

    public static List<PredictionData> loadPredictionData(final String sampleId, final String predictionsDir)
    {
        List<PredictionData> predictionList = Lists.newArrayList();

        try
        {
            String filename = predictionsDir + sampleId + MHCF_PREDICTION_FILE_ID;

            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            // NeId,Allele,Peptide,Affinity=affinity,AffinityPerc=affinity_percentile,
            // ProcScore=processing_score,PresScore=presentation_score,PresPerc=presentation_percentile

            int alleleIndex = fieldsIndexMap.get("HlaAllele");
            int neIdIndex = fieldsIndexMap.get("NeId");
            int peptideIndex = fieldsIndexMap.get("Peptide");
            int affinityIndex = fieldsIndexMap.get("affinity");
            int presentationIndex = fieldsIndexMap.get("presentation_score");

            for(String line : lines)
            {
                predictionList.add(PredictionData.fromCsv(line, alleleIndex, neIdIndex, peptideIndex, affinityIndex, presentationIndex));
            }

            NE_LOGGER.debug("sample({}) loaded {} peptide predictions", sampleId, predictionList.size());
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) predictions file: {}", sampleId, exception.toString());
        }

        return predictionList;
    }

}
