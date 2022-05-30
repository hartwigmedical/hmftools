package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaCommon;

import org.jetbrains.annotations.NotNull;

public final class IsofoxNovelSpliceJunctionLoader {

    private IsofoxNovelSpliceJunctionLoader() {
    }

    @NotNull
    public static List<NovelSpliceJunction> load(@NotNull String isofoxAltSpliceJunctionCsv, @NotNull AltSjCohortData altSjCohortData)
            throws IOException {
        List<NovelSpliceJunction> novelJunctions = Lists.newArrayList();

        BufferedReader fileReader = new BufferedReader(new FileReader(isofoxAltSpliceJunctionCsv));

        String line = fileReader.readLine();
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, RnaCommon.DELIMITER);

        while ((line = fileReader.readLine()) != null) {
            String[] items = line.split(RnaCommon.DELIMITER, -1);

            int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int geneName = fieldsIndexMap.get(FLD_GENE_NAME);
            int chr = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStart = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEnd = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int type = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
            int fragCount = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);
            int depthStart = fieldsIndexMap.get("DepthStart");
            int depthEnd = fieldsIndexMap.get("DepthEnd");
            int regionStart =
                    fieldsIndexMap.containsKey("RegionStart") ? fieldsIndexMap.get("RegionStart") : fieldsIndexMap.get("ContextStart");
            int regionEnd = fieldsIndexMap.containsKey("RegionEnd") ? fieldsIndexMap.get("RegionEnd") : fieldsIndexMap.get("ContextEnd");
            int basesStart = fieldsIndexMap.get("BasesStart");
            int basesEnd = fieldsIndexMap.get("BasesEnd");
            int transStart = fieldsIndexMap.get("TransStart");
            int transEnd = fieldsIndexMap.get("TransEnd");

            AltSpliceJunctionFile altSJ = AltSpliceJunctionFile.fromCsv(items,
                    geneIdIndex,
                    geneName,
                    chr,
                    posStart,
                    posEnd,
                    type,
                    fragCount,
                    depthStart,
                    depthEnd,
                    regionStart,
                    regionEnd,
                    basesStart,
                    basesEnd,
                    transStart,
                    transEnd);

            int cohortFrequency = altSjCohortData.getCohortFrequency(altSJ.key());

            novelJunctions.add(ImmutableNovelSpliceJunction.builder()
                    .geneName(items[fieldsIndexMap.get(FLD_GENE_NAME)])
                    .chromosome(altSJ.Chromosome)
                    .junctionStart(altSJ.SpliceJunction[SE_START])
                    .junctionEnd(altSJ.SpliceJunction[SE_END])
                    .type(altSJ.Type.toString())
                    .fragmentCount(altSJ.FragmentCount)
                    .depthStart(altSJ.DepthCounts[SE_START])
                    .depthEnd(altSJ.DepthCounts[SE_END])
                    .regionStart(altSJ.RegionContexts[SE_START].toString())
                    .regionEnd(altSJ.RegionContexts[SE_END].toString())
                    .basesStart(altSJ.BaseContexts[SE_START])
                    .basesEnd(altSJ.BaseContexts[SE_END])
                    .cohortFrequency(cohortFrequency)
                    .build());
        }
        return novelJunctions;
    }
}
