package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.svs.SvDataType.LINE;
import static com.hartwig.hmftools.cup.svs.SvDataType.MAX_COMPLEX_SIZE;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DEL_20KB_1MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_100KB_5MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.cup.svs.SvDataType.TELOMERIC_SGL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SvDataLoader
{
    public static boolean loadSvDataFromCohortFile(final String filename, final Map<String,SvData> sampleSvData)
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
            CUP_LOGGER.error("failed to read SV cohort data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadSvDataFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,SvData> sampleSvData)
    {
        int i = 0;
        int nextLog = 100;
        for(final String sampleId : sampleIds)
        {
            final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
            final List<LinxCluster> clusterList = dbAccess.readClusters(sampleId);

            mapSvData(sampleId, sampleSvData, svDataList, clusterList);

            ++i;
            if(i >= nextLog)
            {
                nextLog += 100;
                CUP_LOGGER.debug("loaded {} sample SV data sets", i);
            }
        }

        return true;
    }

    public static boolean loadSvDataFromFile(
            final String sampleId, final String svVcfFile, final String clusterFile,
            final Map<String,SvData> sampleSvData)
    {
        try
        {
            final List<StructuralVariantData> svDataList = Lists.newArrayList();
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(svVcfFile, new AlwaysPassFilter());
            final List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;
            for (EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }

            final List<LinxCluster> clusterList = LinxCluster.read(clusterFile);

            mapSvData(sampleId, sampleSvData, svDataList, clusterList);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to load SV file({}) or cluster file({}): {}",
                    sampleId, svVcfFile, clusterFile, e.toString());
            return false;
        }

        return true;
    }

    private static void mapSvData(
            final String sampleId, final Map<String,SvData> sampleSvData,
            final List<StructuralVariantData> allSVs, final List<LinxCluster> clusterList)
    {
        int lineCount = clusterList.stream().filter(x -> x.resolvedType().equals("LINE")).mapToInt(x -> x.clusterCount()).sum();

        // ensure only filtered SVs are considered
        final List<StructuralVariantData> svDataList = allSVs.stream()
                .filter(x -> x.filter().isEmpty() || x.filter().equals(PASS))
                .collect(Collectors.toList());

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

        int maxEventSize = clusterList.stream()
                .filter(x -> !x.resolvedType().equals("LINE"))
                .mapToInt(x -> x.clusterCount()).max().orElse(0);

        SvData svData = new SvData(sampleId);
        svData.setCount(LINE, lineCount);
        svData.setCount(SIMPLE_DEL_20KB_1MB, shortDels);
        svData.setCount(SIMPLE_DUP_32B_200B, shortDups);
        svData.setCount(SIMPLE_DUP_100KB_5MB, longDups);
        svData.setCount(MAX_COMPLEX_SIZE, maxEventSize);
        svData.setCount(TELOMERIC_SGL, telomericSgls);

        sampleSvData.put(sampleId, svData);
    }

    public static boolean loadRefPercentileData(final String filename, final Map<SvDataType,Map<String,double[]>> refSvTypePercentiles)
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
            return false;
        }

        return true;
    }

}
