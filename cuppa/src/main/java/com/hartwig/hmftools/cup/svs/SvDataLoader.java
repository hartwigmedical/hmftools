package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.svs.SvDataType.LINE;
import static com.hartwig.hmftools.cup.svs.SvDataType.MAX_EVENT_SIZE;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DEL_20KB_1MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_100KB_5MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.cup.svs.SvDataType.TELOMERIC_SGL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SvDataLoader
{
    public static void loadSvDataFromCohortFile(final String filename, final Map<String,SvData> sampleSvData)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SvData svData = SvData.from(fieldsIndexMap, line);
                sampleSvData.put(svData.SampleId, svData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read SV data file({}): {}", filename, e.toString());
        }
    }

    public static void loadSvDataFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,SvData> sampleSvData)
    {
        if(dbAccess == null)
            return;

        for(final String sampleId : sampleIds)
        {
            /*'
              summarise(LINE=sum((LEStart!='NONE'&LEStart!='None')|(LEEnd!='NONE'&LEEnd!='None')),
            SIMPLE_DEL_20KB_1MB=sum(Type=='DEL'&Length>=2e4&Length<=1e6),
            SIMPLE_DUP_32B_200B=sum(Type=='DUP'&Length>=32&Length<=200),
            SIMPLE_DUP_100KB_5MB=sum(Type=='DUP'&Length>=1e5&Length<=5e6),
            MAX_EVENT_SIZE=max(ClusterCount),
            TELOMERIC_SGL=sum(Type=='SGL'&(RepeatClass=='Simple_repeat'&(RepeatType=='(TTAGGG)n'|RepeatType=='(CCCTAA)n'))))

            */
            final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);

            // final List<LinxSvAnnotation> svAnnotationList = dbAccess.readSvAnnotaions(sampleId);

            final List<LinxCluster> clusterList = dbAccess.readClusters(sampleId);

            int lineCount = clusterList.stream().filter(x -> x.subType().equals("LINE")).mapToInt(x -> x.clusterCount()).sum();

            int telomericSgls = (int)svDataList.stream()
                    .filter(x -> x.type() == SGL)
                    .filter(x -> x.insertSequenceRepeatClass().equals("Simple_repeat"))
                    .filter(x -> x.insertSequenceRepeatType().equals("(TTAGGG)n") || x.insertSequenceRepeatType().equals("(CCCTAA)n")).count();

            int shortDels = (int)svDataList.stream()
                    .filter(x -> x.type() == DEL)
                    .mapToInt(x -> x.endPosition() - x.startPosition())
                    .filter(x -> x >= 2e4 && x <= 1e6).count();

            int shortDups = (int)svDataList.stream()
                    .filter(x -> x.type() == DUP)
                    .mapToInt(x -> x.endPosition() - x.startPosition())
                    .filter(x -> x >= 32 && x <= 200).count();

            int longDups = (int)svDataList.stream()
                    .filter(x -> x.type() == DUP)
                    .mapToInt(x -> x.endPosition() - x.startPosition())
                    .filter(x -> x >= 1e5 && x <= 5e6).count();

            int maxEventSize = clusterList.stream().mapToInt(x -> x.clusterCount()).max().orElse(0);

            SvData svData = new SvData(sampleId);
            svData.setCount(LINE, lineCount);
            svData.setCount(SIMPLE_DEL_20KB_1MB, shortDels);
            svData.setCount(SIMPLE_DUP_32B_200B, shortDups);
            svData.setCount(SIMPLE_DUP_100KB_5MB, longDups);
            svData.setCount(MAX_EVENT_SIZE, maxEventSize);
            svData.setCount(TELOMERIC_SGL, telomericSgls);

            sampleSvData.put(sampleId, svData);
        }
    }

    public static void loadRefPercentileData(final String filename, final Map<SvDataType,Map<String,double[]>> refSvTypePercentiles)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String cancerType = items[0];
                final SvDataType svDataType = SvDataType.valueOf(items[1]);

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                Map<String,double[]> svPercData = refSvTypePercentiles.get(svDataType);

                if(svPercData == null)
                {
                    svPercData = Maps.newHashMap();
                    refSvTypePercentiles.put(svDataType, svPercData);
                }

                svPercData.put(cancerType, percentileData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read SV perc data file({}): {}", filename, e.toString());
        }
    }

}
