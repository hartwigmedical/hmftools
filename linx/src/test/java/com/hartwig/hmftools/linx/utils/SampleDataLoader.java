package com.hartwig.hmftools.linx.utils;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.utils.SvTestUtils.initialiseSV;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFile;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.NotNull;

public class SampleDataLoader
{
    // data input expects all StructuralVariant data fields, followed by additional svAnnotation fields

    private static final int COL_PLOIDY_MIN = 53;
    private static final int COL_PLOIDY_MAX = 54;
    private static final int COL_EXPECTED = COL_PLOIDY_MAX + 1;

    public static List<SvVarData> loadSampleTestData(@NotNull final String resource)
    {
        final InputStream inputStream = StructuralVariantFile.class.getResourceAsStream("/sample_data/" + resource);

        final List<String> inputData = new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("svId"))
                .collect(toList());

        final List<StructuralVariantData> svData = inputData.stream().map(StructuralVariantFile::fromString).collect(toList());

        return createSVs(svData, inputData);
    }

    private static List<SvVarData> createSVs(final List<StructuralVariantData> svDataList, final List<String> inputData)
    {
        List<SvVarData> svList = Lists.newArrayList();

        for(int i = 0; i < svDataList.size(); ++i)
        {
            final StructuralVariantData svData = svDataList.get(i);
            SvVarData var = new SvVarData(svData);
            svList.add(var);
            initialiseSV(var);

            final String[] items = inputData.get(i).split(StructuralVariantFile.DELIMITER);

            if(items.length < COL_EXPECTED)
                return svList;

            var.setJcnRecalcData(Double.parseDouble(items[COL_PLOIDY_MIN]), Double.parseDouble(items[COL_PLOIDY_MAX]));
        }

        return svList;
    }

}
