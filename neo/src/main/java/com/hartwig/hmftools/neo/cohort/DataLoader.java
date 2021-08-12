package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.extractTranscriptNames;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.IM_FILE_ID;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.BindingPredictionData.DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.utils.Matrix;

public class DataLoader
{
    public static final String LILAC_COVERAGE_FILE_ID = ".lilac.csv";
    public static final String MCF_PREDICTION_FILE_ID = ".mcf.predictions.csv";
    public static final String MCF_AFFINITY = "affinity";
    public static final String MCF_AFFINITY_PERC = "affinity_percentile";
    public static final String MCF_PRESENTATION = "presentation_score";
    public static final String MCF_PRESENTATION_PERC = "presentation_percentile";

    public static Map<Integer,NeoEpitopeData> loadNeoEpitopes(final String sampleId, final String neoDataDir)
    {
        Map<Integer,NeoEpitopeData> neoDataMap = Maps.newHashMap();

        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(neoDataDir, sampleId);
            String rnaNeoEpitopeFile = neoEpitopeFile.replace(IM_FILE_ID, ISF_FILE_ID);

            boolean hasRnaData = Files.exists(Paths.get(rnaNeoEpitopeFile));

            final List<String> lines = hasRnaData ?
                    Files.readAllLines(new File(rnaNeoEpitopeFile).toPath()) : Files.readAllLines(new File(neoEpitopeFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int neIdIndex = fieldsIndexMap.get("NeId");
            int varTypeIndex = fieldsIndexMap.get("VariantType");
            int varInfoIndex = fieldsIndexMap.get("VariantInfo");
            int geneNameUpIndex = fieldsIndexMap.get("GeneNameUp");
            int geneNameDownIndex = fieldsIndexMap.get("GeneNameDown");
            int geneIdDownIndex = fieldsIndexMap.get("GeneIdDown");
            int tpmCancerIndex = fieldsIndexMap.get("TpmCancerDown");
            int tpmCohortIndex = fieldsIndexMap.get("TpmCohortDown");
            int upAaIndex = fieldsIndexMap.get("UpstreamAA");
            int downAaIndex = fieldsIndexMap.get("DownstreamAA");
            int novelAaIndex = fieldsIndexMap.get("NovelAA");
            int transUpIndex = fieldsIndexMap.get("UpTranscripts");
            int transDownIndex = fieldsIndexMap.get("DownTranscripts");

            Integer rnaFragIndex = fieldsIndexMap.get("FragmentsNovel");
            Integer rnaBaseDepthUp = fieldsIndexMap.get("BaseDepthUp");
            Integer rnaBaseDepthDown = fieldsIndexMap.get("BaseDepthDown");

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                int neId = Integer.parseInt(items[neIdIndex]);

                String geneNameUp = items[geneNameUpIndex];
                String geneNameDown = items[geneNameDownIndex];
                String geneIdDown = items[geneIdDownIndex];
                String geneName = geneNameUp.equals(geneNameDown) ? geneNameUp : geneNameUp + "_" + geneNameDown;

                String aminoAcids = items[upAaIndex] + items[novelAaIndex] + items[downAaIndex];

                double tpmCancer = Double.parseDouble(items[tpmCancerIndex]);
                double tpmCohort = Double.parseDouble(items[tpmCohortIndex]);

                List<String> transUpNames = Lists.newArrayList();
                List<String> transDownNames = Lists.newArrayList();

                extractTranscriptNames(items[transUpIndex], items[transDownIndex], transUpNames, transDownNames);

                int rnaFragCount = 0;
                int[] rnaBaseDepth = new int[] {0, 0};

                if(hasRnaData)
                {
                    rnaFragCount = Integer.parseInt(items[rnaFragIndex]);
                    rnaBaseDepth[SE_START] =  Integer.parseInt(items[rnaBaseDepthUp]);
                    rnaBaseDepth[SE_END] =  Integer.parseInt(items[rnaBaseDepthDown]);
                }

                NeoEpitopeData neoData = new NeoEpitopeData(
                        neId, NeoEpitopeType.valueOf(items[varTypeIndex]), items[varInfoIndex], geneIdDown, geneName,
                        aminoAcids, transUpNames, transDownNames, tpmCancer, tpmCohort, rnaFragCount, rnaBaseDepth);

                neoDataMap.put(neId, neoData);
            }

            NE_LOGGER.debug("sample({}) loaded {} neo-epitopes", sampleId, lines.size());
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) neo-epitope file: {}", sampleId, exception.toString());
        }

        return neoDataMap;
    }


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

    public static List<BindingPredictionData> loadPredictionData(final String sampleId, final String predictionsDir)
    {
        List<BindingPredictionData> predictionList = Lists.newArrayList();

        try
        {
            String filename = predictionsDir + sampleId + MCF_PREDICTION_FILE_ID;

            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            // NeId,Allele,Peptide,Affinity=affinity,AffinityPerc=affinity_percentile,
            // ProcScore=processing_score,PresScore=presentation_score,PresPerc=presentation_percentile

            int alleleIndex = fieldsIndexMap.get("HlaAllele");
            int neIdIndex = fieldsIndexMap.get("NeId");
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            int affinityIndex = fieldsIndexMap.get(MCF_AFFINITY);
            int affinityPercIndex = fieldsIndexMap.get(MCF_AFFINITY_PERC);
            int presIndex = fieldsIndexMap.get(MCF_PRESENTATION);
            int presPercIndex = fieldsIndexMap.get(MCF_PRESENTATION_PERC);

            for(String line : lines)
            {
                predictionList.add(BindingPredictionData.fromMcfCsv(
                        line, alleleIndex, neIdIndex, peptideIndex, affinityIndex, affinityPercIndex, presIndex, presPercIndex));
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
