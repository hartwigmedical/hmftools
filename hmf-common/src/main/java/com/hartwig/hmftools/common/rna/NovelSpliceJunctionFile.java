package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_END;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_DEPTH_START;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public final class NovelSpliceJunctionFile
{
    public static final String ALT_SJ_FILE_ID = "alt_splice_junc.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + ALT_SJ_FILE_ID;
    }

    public static final String FLD_ALT_SJ_POS_START = "SjStart";
    public static final String FLD_ALT_SJ_POS_END = "SjEnd";
    public static final String FLD_ALT_SJ_TYPE = "Type";
    public static final String FLD_BASES_START = "BaseStart";
    public static final String FLD_BASES_END = "BaseEnd";
    public static final String FLD_COHORT_FREQUENCY = "CohortFrequency";

    public static List<NovelSpliceJunction> read(final String filename)
    {
        try
        {
            List<NovelSpliceJunction> novelJunctions = Lists.newArrayList();

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            String fileDelim = inferFileDelimiter(filename);
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, fileDelim);

            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim, -1);

                int geneIndex = fieldsIndexMap.get(FLD_GENE_NAME);
                int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
                int sjStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
                int sjEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
                int typeIndex = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
                int fragCountIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);
                int depthStartIndex = fieldsIndexMap.get(FLD_DEPTH_START);
                int depthEndIndex = fieldsIndexMap.get(FLD_DEPTH_END);
                int regionStartIndex = fieldsIndexMap.get(FLD_REGION_START);
                int regionEndIndex = fieldsIndexMap.get(FLD_REGION_END);
                Integer basesStartIndex = fieldsIndexMap.get(FLD_BASES_START);
                Integer basesEndIndex = fieldsIndexMap.get(FLD_BASES_END);
                Integer cohortFreqIndex = fieldsIndexMap.get(FLD_COHORT_FREQUENCY);

                // int transStart = fieldsIndexMap.get("TransStart");
                // int transEnd = fieldsIndexMap.get("TransEnd");

                String chromosome = values[chrIndex];
                int sjStart = Integer.parseInt(values[sjStartIndex]);
                int sjEnd = Integer.parseInt(values[sjEndIndex]);


                novelJunctions.add(ImmutableNovelSpliceJunction.builder()
                        .geneName(values[geneIndex])
                        .chromosome(chromosome)
                        .junctionStart(sjStart)
                        .junctionEnd(sjEnd)
                        .type(AltSpliceJunctionType.valueOf(values[typeIndex]))
                        .fragmentCount(Integer.parseInt(values[fragCountIndex]))
                        .depthStart(Integer.parseInt(values[depthStartIndex]))
                        .depthEnd(Integer.parseInt(values[depthEndIndex]))
                        .regionStart(AltSpliceJunctionContext.valueOf(values[regionStartIndex]))
                        .regionEnd(AltSpliceJunctionContext.valueOf(values[regionEndIndex]))
                        .basesStart(basesStartIndex != null ? values[basesStartIndex] : "")
                        .basesEnd(basesEndIndex != null ? values[basesEndIndex] : "")
                        .cohortFrequency(cohortFreqIndex != null ? Integer.parseInt(values[cohortFreqIndex]) : 0)
                        .build());
            }

            return novelJunctions;
        }
        catch(IOException e)
        {
            RNA_LOGGER.error("failed to load Isofox fusion file({}): {}", filename, e.toString());
            return null;
        }
    }
}
