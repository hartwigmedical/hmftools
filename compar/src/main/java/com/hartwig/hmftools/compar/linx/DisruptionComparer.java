package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_CODING_CONTEXT;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_GENE_ORIENT;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_NEXT_SPLICE;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_REGION_TYPE;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DisruptionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DisruptionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
    }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_REGION_TYPE, FLD_CODING_CONTEXT, FLD_GENE_ORIENT, FLD_NEXT_SPLICE);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
        final List<LinxBreakend> breakends = dbAccess.readBreakends(sampleId);
        return buildBreakends(svDataList, breakends);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            final List<StructuralVariantData> svDataList = Lists.newArrayList();

            String vcfFile = PurpleCommon.purpleSvFile(fileSources.Purple, sampleId);
            List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());
            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;

            for(EnrichedStructuralVariant variant : enrichedVariants)
            {
                svDataList.add(convertSvData(variant, svId++)); // valid to set ID again since read this way in Linx
            }

            List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(fileSources.Linx, sampleId));

            CMP_LOGGER.debug("sample({}) loaded {} SVs {} breakends",sampleId, svDataList.size(), breakends.size());

            return buildBreakends(svDataList, breakends);

        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx breakend disruption data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }

    private List<ComparableItem> buildBreakends(final List<StructuralVariantData> svDataList, final List<LinxBreakend> breakends)
    {
        List<ComparableItem> items = Lists.newArrayList();

        for(StructuralVariantData var : svDataList)
        {
            List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == var.id()).collect(Collectors.toList());

            for(LinxBreakend breakend : svBreakends)
            {
                breakends.remove(breakend);

                DisruptionData disruptionData = new DisruptionData(var, breakend);
                items.add(disruptionData);
            }
        }

        return items;

    }
}
