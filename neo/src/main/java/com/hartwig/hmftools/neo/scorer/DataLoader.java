package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_AA_DOWN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_AA_NOVEL;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_AA_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_CB_LEN_MAX;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_CB_LEN_MIN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_FUSED_LEN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_GENE_ID_DOWN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_GENE_NAME_DOWN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_GENE_NAME_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_ID;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_CN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_NMD_MAX;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_NMD_MIN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_SC_LIKELIHOOD;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_SKIP_ACCEPTORS;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_SKIP_DONORS;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_TRANS_DOWN;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_TRANS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_VAR_INFO;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_VAR_TYPE;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.extractTranscriptNames;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.scorer.AlleleCoverage.EXPECTED_ALLELE_COUNT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class DataLoader
{
    public static final String LILAC_COVERAGE_FILE_ID = ".lilac.csv";

    public static List<NeoEpitopeData> loadNeoEpitopes(final String sampleId, final String neoDataDir)
    {
        // load into lists of neoepitopes by variant
        List<NeoEpitopeData> neoDataList = Lists.newArrayList();

        try
        {
            String neoEpitopeFile = NeoEpitopeFile.generateFilename(neoDataDir, sampleId);

            if(!Files.exists(Paths.get(neoEpitopeFile)))
            {
                neoEpitopeFile = neoDataDir + sampleId + ".imu.neo_epitopes.csv";
            }

            final List<String> lines = Files.readAllLines(new File(neoEpitopeFile).toPath());

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int neIdIndex = fieldsIndexMap.get(FLD_NE_ID);
            int varTypeIndex = fieldsIndexMap.get(FLD_NE_VAR_TYPE);
            int varInfoIndex = fieldsIndexMap.get(FLD_NE_VAR_INFO);
            int geneIdUpIndex = fieldsIndexMap.get(FLD_NE_GENE_NAME_UP);
            int geneIdDownIndex = fieldsIndexMap.get(FLD_NE_GENE_ID_DOWN);
            int geneNameUpIndex = fieldsIndexMap.get(FLD_NE_GENE_NAME_UP);
            int geneNameDownIndex = fieldsIndexMap.get(FLD_NE_GENE_NAME_DOWN);
            int upAaIndex = fieldsIndexMap.get(FLD_NE_AA_UP);
            int downAaIndex = fieldsIndexMap.get(FLD_NE_AA_DOWN);
            int novelAaIndex = fieldsIndexMap.get(FLD_NE_AA_NOVEL);
            int nmdMinIndex = fieldsIndexMap.get(FLD_NE_NMD_MIN);
            int nmdMaxIndex = fieldsIndexMap.get(FLD_NE_NMD_MAX);
            int cnIndex = fieldsIndexMap.get(FLD_NE_CN);
            int sclIndex = fieldsIndexMap.get(FLD_NE_SC_LIKELIHOOD);
            int cbLenMinIndex = fieldsIndexMap.get(FLD_NE_CB_LEN_MIN);
            int cbLenMaxIndex = fieldsIndexMap.get(FLD_NE_CB_LEN_MAX);
            int feLenIndex = fieldsIndexMap.get(FLD_NE_FUSED_LEN);
            int skipDonIndex = fieldsIndexMap.get(FLD_NE_SKIP_DONORS);
            int skipAccIndex = fieldsIndexMap.get(FLD_NE_SKIP_ACCEPTORS);
            int transUpIndex = fieldsIndexMap.get(FLD_NE_TRANS_UP);
            int transDownIndex = fieldsIndexMap.get(FLD_NE_TRANS_DOWN);

            for(String line : lines)
            {
                final String[] values = line.split(DELIMITER, -1);

                int neId = Integer.parseInt(values[neIdIndex]);

                String geneNameUp = values[geneNameUpIndex];
                String geneNameDown = values[geneNameDownIndex];
                String geneIdDown = values[geneIdDownIndex];
                String geneName = geneNameUp.equals(geneNameDown) ? geneNameUp : geneNameUp + "_" + geneNameDown;

                List<String> transUpNames = Lists.newArrayList();
                List<String> transDownNames = Lists.newArrayList();

                extractTranscriptNames(values[transUpIndex], values[transDownIndex], transUpNames, transDownNames);

                NeoEpitopeData neoData = new NeoEpitopeData(
                        neId, NeoEpitopeType.valueOf(values[varTypeIndex]), values[varInfoIndex], geneIdDown, geneName,
                        values[upAaIndex],values[novelAaIndex], values[downAaIndex], transUpNames, transDownNames,
                        Integer.parseInt(values[nmdMinIndex]), Integer.parseInt(values[nmdMaxIndex]),
                        Double.parseDouble(values[cnIndex]), Double.parseDouble(values[sclIndex]),
                        Integer.parseInt(values[cbLenMinIndex]), Integer.parseInt(values[cbLenMaxIndex]),
                        Integer.parseInt(values[feLenIndex]), Integer.parseInt(values[skipDonIndex]), Integer.parseInt(values[skipAccIndex]));

                neoDataList.add(neoData);
            }

            NE_LOGGER.debug("sample({}) loaded {} neo-epitopes", sampleId, neoDataList.size());
        }
        catch(IOException exception)
        {
            NE_LOGGER.error("failed to read sample({}) neo-epitope file: {}", sampleId, exception.toString());
        }

        return neoDataList;
    }

    public static List<SomaticVariant> loadSomaticVariants(
            final String sampleId, final String rnaSampleId, final String rnaSomaticVcf, final List<NeoEpitopeData> pointNeos)
    {
        List<SomaticVariant> matchedVariants = Lists.newArrayList();

        String vcfFile = rnaSomaticVcf.replaceAll("\\*", sampleId);

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;

                SomaticVariant variant = SomaticVariant.fromContext(variantContext, sampleId, rnaSampleId);

                if(variant == null)
                    continue;

                String varInfo = variant.variantInfo();

                if(pointNeos.stream().anyMatch(x -> x.VariantInfo.equals(varInfo)))
                {
                    matchedVariants.add(variant);
                }
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return matchedVariants;
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

            List<RnaNeoEpitope> rnaNeoDataList = RnaNeoEpitope.read(filename).stream()
                    .filter(x -> x.VariantType.isFusion()).collect(Collectors.toList());

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
