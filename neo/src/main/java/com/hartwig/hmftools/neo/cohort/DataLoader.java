package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.extractTranscriptNames;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;

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
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;

public class DataLoader
{
    public static final String LILAC_COVERAGE_FILE_ID = ".lilac.csv";

    public static Map<Integer,NeoEpitopeData> loadNeoEpitopes(final String sampleId, final String neoDataDir)
    {
        Map<Integer,NeoEpitopeData> neoDataMap = Maps.newHashMap();

        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(neoDataDir, sampleId);

            final List<String> lines = Files.readAllLines(new File(neoEpitopeFile).toPath());

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

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                int neId = Integer.parseInt(items[neIdIndex]);

                String geneNameUp = items[geneNameUpIndex];
                String geneNameDown = items[geneNameDownIndex];
                String geneIdDown = items[geneIdDownIndex];
                String geneName = geneNameUp.equals(geneNameDown) ? geneNameUp : geneNameUp + "_" + geneNameDown;

                double tpmCancer = Double.parseDouble(items[tpmCancerIndex]);
                double tpmCohort = Double.parseDouble(items[tpmCohortIndex]);

                List<String> transUpNames = Lists.newArrayList();
                List<String> transDownNames = Lists.newArrayList();

                extractTranscriptNames(items[transUpIndex], items[transDownIndex], transUpNames, transDownNames);

                NeoEpitopeData neoData = new NeoEpitopeData(
                        neId, NeoEpitopeType.valueOf(items[varTypeIndex]), items[varInfoIndex], geneIdDown, geneName,
                        items[upAaIndex],items[novelAaIndex], items[downAaIndex], transUpNames, transDownNames, tpmCancer, tpmCohort);

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

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
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

    public static List<RnaNeoEpitope> loadRnaNeoData(final String sampleId, final String isofoxDataDir)
    {
        if(isofoxDataDir == null)
            return Lists.newArrayList();

        try
        {
            String filename = RnaNeoEpitope.generateFilename(isofoxDataDir, sampleId);

            if(!Files.exists(Paths.get(filename)))
                return Lists.newArrayList();

            List<RnaNeoEpitope> rnaNeoDataList = RnaNeoEpitope.read(filename);

            NE_LOGGER.debug("sample({}) loaded {} RNA neoepitopes", sampleId, rnaNeoDataList.size());

            return rnaNeoDataList;
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) predictions file: {}", sampleId, exception.toString());
            return Lists.newArrayList();
        }
    }

}
