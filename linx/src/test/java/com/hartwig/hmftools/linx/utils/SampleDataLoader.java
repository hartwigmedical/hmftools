package com.hartwig.hmftools.linx.utils;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.initialiseSV;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;
import com.hartwig.hmftools.linx.gene.BreakendGenePrep;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.NotNull;

public class SampleDataLoader
{
    // data input expects all StructuralVariant data fields, followed by additional svAnnotation fields

    private static final String PLOIDY_MIN = "ploidyMin";
    private static final String PLOIDY_MAX = "ploidyMax";

    private static final int COL_PLOIDY_MIN = 53;
    private static final int COL_PLOIDY_MAX = 54;
    private static final int COL_EXPECTED = COL_PLOIDY_MAX + 1;

    public static List<SvVarData> loadSampleTestData(@NotNull final String resource)
    {
        InputStream inputStream = SampleDataLoader.class.getResourceAsStream("/sample_data/" + resource);

        List<String> lines = new BufferedReader(new InputStreamReader(inputStream)).lines()
                .collect(toList());

        List<SvVarData> svList = Lists.newArrayList();

        String header = lines.get(0);
        lines.remove(0);
        Map<String,Integer> fieldsIndexMap = FileReaderUtils.createFieldsIndexMap(header, TSV_DELIM);

        int ploidyMinIndex = fieldsIndexMap.get(PLOIDY_MIN);
        int ploidyMaxIndex = fieldsIndexMap.get(PLOIDY_MAX);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            StructuralVariantData svData = StructuralVariantFile.fromString(values);

            SvVarData var = new SvVarData(svData);
            svList.add(var);
            initialiseSV(var);

            double ploidyMin = Double.parseDouble(values[ploidyMinIndex]);
            double ploidyMax = Double.parseDouble(values[ploidyMaxIndex]);
            var.setJcnRecalcData(ploidyMin, ploidyMax);
        }

        return svList;
    }

    private static List<SvVarData> createSVs(final List<StructuralVariantData> svDataList, final List<String> inputData)
    {
        List<SvVarData> svList = Lists.newArrayList();

        for(int i = 0; i < svDataList.size(); ++i)
        {
            StructuralVariantData svData = svDataList.get(i);
            SvVarData var = new SvVarData(svData);
            svList.add(var);
            initialiseSV(var);

            final String[] items = inputData.get(i).split(TSV_DELIM);

            if(items.length < COL_EXPECTED)
                return svList;

            var.setJcnRecalcData(Double.parseDouble(items[COL_PLOIDY_MIN]), Double.parseDouble(items[COL_PLOIDY_MAX]));
        }

        return svList;
    }

    public static void setSvGeneData(final List<SvVarData> svList, final EnsemblDataCache ensemblDataCache, boolean applyPreGeneDistance)
    {
        int preGeneDistance = applyPreGeneDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;
        BreakendGenePrep.setSvGeneData(svList, ensemblDataCache, preGeneDistance, Maps.newHashMap());
    }
}
