package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.cup.prep.CategoryType.SV;
import static com.hartwig.hmftools.cup.svs.SvDataType.LINE;
import static com.hartwig.hmftools.cup.svs.SvDataType.MAX_COMPLEX_SIZE;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DEL_20KB_1MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_100KB_5MB;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.cup.svs.SvDataType.TELOMERIC_SGL;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;

import java.io.File;
import java.nio.file.NoSuchFileException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class StructuralVariantPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    public StructuralVariantPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return SV; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        List<DataItem> dataItems = Lists.newArrayList();

        try
        {
            final List<StructuralVariantData> svDataList = Lists.newArrayList();

            String purpleSvFile = mConfig.purpleSvFile(sampleId);
            if(!new File(purpleSvFile).isFile())
                throw new NoSuchFileException(purpleSvFile);
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(purpleSvFile, new AlwaysPassFilter());
            final List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;
            for (EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }

            String linxClusterFile = mConfig.linxClusterFile(sampleId);
            if(!new File(linxClusterFile).isFile())
                throw new NoSuchFileException(linxClusterFile);
            final List<LinxCluster> clusterList = LinxCluster.read(linxClusterFile);

            final int[] svDataCounts = extractSvCounts(svDataList, clusterList);

            for(SvDataType type : SvDataType.values())
            {
                dataItems.add(new DataItem(DNA, ItemType.SV_COUNT, type.toString(), svDataCounts[type.ordinal()]));
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to extract category({}):", sampleId, categoryType());
            e.printStackTrace();
            System.exit(1);
        }

        return dataItems;
    }

    public static int[] extractSvCounts(final List<StructuralVariantData> allSVs, final List<LinxCluster> clusterList)
    {
        final int[] svCounts = new int[SvDataType.values().length];

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

        svCounts[LINE.ordinal()] = lineCount;
        svCounts[SIMPLE_DEL_20KB_1MB.ordinal()] = shortDels;
        svCounts[SIMPLE_DUP_32B_200B.ordinal()] = shortDups;
        svCounts[SIMPLE_DUP_100KB_5MB.ordinal()] = longDups;
        svCounts[MAX_COMPLEX_SIZE.ordinal()] = maxEventSize;
        svCounts[TELOMERIC_SGL.ordinal()] = telomericSgls;

        return svCounts;
    }
}
