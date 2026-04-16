package com.hartwig.hmftools.compar.linx;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.compar.common.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_BREAKEND_INFO;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import htsjdk.tribble.TribbleException;

public class DisruptionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    private final Map<SourceType,List<LinxBreakend>> mBreakends;
    private final Map<SourceType,List<StructuralVariantData>> mSvDataList;

    public DisruptionComparer(final ComparConfig config)
    {
        mConfig = config;
        mBreakends = Maps.newHashMap();
        mSvDataList = Maps.newHashMap();
    }

    public Map<SourceType,List<LinxBreakend>> breakends() { return mBreakends; }
    public Map<SourceType,List<StructuralVariantData>> svDataList() { return mSvDataList; }

    @Override
    public CategoryType category() { return DISRUPTION; }

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
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        List<StructuralVariantData> svDataList = Lists.newArrayList(dbAccess.readStructuralVariantData(sampleId));
        List<LinxBreakend> breakends = Lists.newArrayList(dbAccess.readBreakends(sampleId));

        mSvDataList.put(sourceType, svDataList);
        mBreakends.put(sourceType, breakends);

        return buildBreakends(sourceType);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<StructuralVariantData> svDataList = Lists.newArrayList();
        List<LinxBreakend> breakends = Lists.newArrayList();

        mSvDataList.put(fileSources.Source, svDataList);
        mBreakends.put(fileSources.Source, breakends);

        try
        {
            String vcfFile = PurpleCommon.purpleSomaticSvFile(fileSources.Purple, sampleId);

            List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());

            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;

            for(EnrichedStructuralVariant variant : enrichedVariants)
            {
                svDataList.add(convertSvData(variant, svId++)); // valid to set ID again since read this way in Linx
            }

            breakends.addAll(LinxBreakend.read(LinxBreakend.generateFilename(fileSources.Linx, sampleId)));

            CMP_LOGGER.debug("sample({}) loaded {} SVs {} breakends",sampleId, mSvDataList.size(), mBreakends.size());

            return buildBreakends(fileSources.Source);

        }
        catch(IOException | TribbleException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx breakend disruption data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> buildBreakends(final SourceType sourceType)
    {
        List<ComparableItem> items = Lists.newArrayList();

        Map<String,List<BreakendData>> geneBreakendMap = Maps.newHashMap();

        MatchLevel matchLevel = mConfig.Categories.getOrDefault(DISRUPTION, MatchLevel.REPORTABLE);

        List<StructuralVariantData> svDataList = mSvDataList.get(sourceType);
        List<LinxBreakend> breakends = mBreakends.get(sourceType);

        for(StructuralVariantData var : svDataList)
        {
            List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == var.id()).collect(Collectors.toList());

            for(LinxBreakend breakend : svBreakends)
            {
                // breakends.remove(breakend); was an optimisation

                if(matchLevel == MatchLevel.REPORTABLE && breakend.reportedStatus() != ReportedStatus.REPORTED)
                    continue;

                List<BreakendData> geneBreakends = geneBreakendMap.get(breakend.gene());

                if(geneBreakends == null)
                {
                    geneBreakends = Lists.newArrayList();
                    geneBreakendMap.put(breakend.gene(), geneBreakends);
                }

                boolean usesStart = breakend.isStart();

                String chromosome = usesStart ? var.startChromosome() : var.endChromosome();
                int position = usesStart ? var.startPosition() : var.endPosition();

                BasePosition comparisonPosition = determineComparisonGenomePosition(
                        chromosome, position, sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                BreakendData breakendData = buildBreakendData(breakend, var, comparisonPosition);

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

    protected static BreakendData buildBreakendData(
            final LinxBreakend breakend, final StructuralVariantData var, @Nullable final BasePosition comparisonPosition)
    {
        boolean usesStart = breakend.isStart();

        int[] homologyOffsets = usesStart ?
                new int[] { var.startIntervalOffsetStart(), var.startIntervalOffsetEnd() } :
                new int[] { var.endIntervalOffsetStart(), var.endIntervalOffsetEnd() };

        String chromosome = usesStart ? var.startChromosome() : var.endChromosome();
        int position = usesStart ? var.startPosition() : var.endPosition();

        int depthStart = var.startTumorReferenceFragmentCount() + var.startNormalReferenceFragmentCount();
        int fragsStart = var.startTumorReferenceFragmentCount();
        int depthEnd = var.endTumorReferenceFragmentCount() + var.endNormalReferenceFragmentCount();
        int fragsEnd = var.endTumorReferenceFragmentCount();
        int qual = (int)round(var.qualityScore());

        return new BreakendData(
            breakend, usesStart ? var.vcfIdStart() : var.vcfIdEnd(), var.type(), chromosome, position,
            usesStart ? var.startOrientation() : var.endOrientation(), homologyOffsets,
            usesStart ? depthStart : depthEnd, usesStart ? fragsStart : fragsEnd, qual,
                comparisonPosition != null ? comparisonPosition.Chromosome : chromosome,
                comparisonPosition != null ? comparisonPosition.Position : position);
    }
}
