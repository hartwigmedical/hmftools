package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantFileLoader.fromGridssFile;
import static com.hartwig.hmftools.common.sv.SvFactoryInterface.buildSvFactory;
import static com.hartwig.hmftools.compar.common.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_BREAKEND_INFO;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.SvFactoryInterface;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import htsjdk.variant.variantcontext.filter.VariantContextFilter;

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
    public void registerThresholds(final DiffThresholds thresholds) {}

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_REPORTED, FLD_BREAKEND_INFO);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
        final List<LinxBreakend> breakends = dbAccess.readBreakends(sampleId);
        return buildBreakends(svDataList, breakends, sourceName);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        try
        {
            final List<StructuralVariantData> svDataList = Lists.newArrayList();

            String vcfFile = PurpleCommon.purpleSomaticSvFile(fileSources.Purple, sampleId);

            List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());

            if(variants.stream().anyMatch(x -> x.startContext().getID().contains("gridss")))
            {
                variants = fromGridssFile(vcfFile, new AlwaysPassFilter());
            }

            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;

            for(EnrichedStructuralVariant variant : enrichedVariants)
            {
                svDataList.add(convertSvData(variant, svId++)); // valid to set ID again since read this way in Linx
            }

            List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(fileSources.Linx, sampleId));

            CMP_LOGGER.debug("sample({}) loaded {} SVs {} breakends",sampleId, svDataList.size(), breakends.size());

            return buildBreakends(svDataList, breakends, fileSources.Source);

        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx breakend disruption data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> buildBreakends(final List<StructuralVariantData> svDataList, final List<LinxBreakend> breakends, final String sourceName)
    {
        List<ComparableItem> items = Lists.newArrayList();

        Map<String,List<BreakendData>> geneBreakendMap = Maps.newHashMap();

        MatchLevel matchLevel = mConfig.Categories.getOrDefault(DISRUPTION, MatchLevel.REPORTABLE);

        for(StructuralVariantData var : svDataList)
        {
            List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == var.id()).collect(Collectors.toList());

            for(LinxBreakend breakend : svBreakends)
            {
                breakends.remove(breakend);

                if(matchLevel == MatchLevel.REPORTABLE && !breakend.reportedDisruption())
                    continue;

                List<BreakendData> geneBreakends = geneBreakendMap.get(breakend.gene());

                if(geneBreakends == null)
                {
                    geneBreakends = Lists.newArrayList();
                    geneBreakendMap.put(breakend.gene(), geneBreakends);
                }

                boolean usesStart = breakend.isStart();

                int[] homologyOffsets = usesStart ?
                        new int[] { var.startIntervalOffsetStart(), var.startIntervalOffsetEnd() } :
                        new int[] { var.endIntervalOffsetStart(), var.endIntervalOffsetEnd() };

                String chromosome = usesStart ? var.startChromosome() : var.endChromosome();
                int position = usesStart ? var.startPosition() : var.endPosition();

                BasePosition comparisonPosition = determineComparisonGenomePosition(
                        chromosome, position, sourceName, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                BreakendData breakendData = new BreakendData(
                        breakend, usesStart ? var.vcfIdStart() : var.vcfIdEnd(), var.type(), chromosome, position,
                        usesStart ? var.startOrientation() : var.endOrientation(), homologyOffsets,
                        comparisonPosition.Chromosome, comparisonPosition.Position);

                geneBreakends.add(breakendData);
            }
        }

        for(Map.Entry<String,List<BreakendData>> entry : geneBreakendMap.entrySet())
        {
            String geneName = entry.getKey();
            List<BreakendData> geneBreakends = entry.getValue();

            DisruptionData disruptionData = new DisruptionData(DISRUPTION, geneName, geneBreakends);
            items.add(disruptionData);
        }

        return items;

    }
}
