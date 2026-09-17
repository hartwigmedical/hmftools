package com.hartwig.hmftools.compar.linx;

import static java.lang.Math.round;

import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_SV;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GermlineSvComparer extends ItemComparer
{
    protected enum Fields
    {
        Reported,
        BreakendInfo;
    }

    public GermlineSvComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Reported.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Reported.toString()), null));

        mFields.add(new FieldInfo(
                Fields.BreakendInfo.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.BreakendInfo.toString()), null));
    }

    @Override
    public CategoryType category() { return GERMLINE_SV; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        List<ComparableItem> items = Lists.newArrayList();

        Map<String,List<BreakendData>> geneBreakendMap = Maps.newHashMap();

        MatchLevel matchLevel = mConfig.MatchingLevel;

        try
        {
            String germlineSvFile = LinxGermlineDisruption.generateFilename(fileSources.LinxGermline, sampleId);
            List<LinxGermlineDisruption> germlineSvs = null; // loads on demand

            String germlineBreakendFile = LinxBreakend.generateFilename(fileSources.LinxGermline, sampleId, true);

            boolean reportedOnly = matchLevel == MatchLevel.REPORTABLE;

            for(LinxBreakend breakend : LinxBreakend.read(germlineBreakendFile))
            {
                if(reportedOnly && breakend.reportedStatus() != ReportedStatus.REPORTED)
                    continue;

                if(germlineSvs == null)
                {
                    germlineSvs = LinxGermlineDisruption.read(germlineSvFile);
                    CMP_LOGGER.debug("sample({}) loaded {} germline SVs", sampleId, germlineSvs.size());
                }

                LinxGermlineDisruption var = germlineSvs.stream().filter(x -> x.SvId == breakend.svId()).findFirst().orElse(null);

                if(var == null)
                    continue; // implies an inconsistency

                List<BreakendData> geneBreakends = geneBreakendMap.get(breakend.gene());

                if(geneBreakends == null)
                {
                    geneBreakends = Lists.newArrayList();
                    geneBreakendMap.put(breakend.gene(), geneBreakends);
                }

                boolean usesStart = breakend.isStart();

                int[] homologyOffsets = {0, 0}; // not available at the moment unless VCF is read

                String chromosome = usesStart ? var.ChromosomeStart : var.ChromosomeEnd;
                int position = usesStart ? var.PositionStart : var.PositionEnd;

                if(mConfig.RequiresLiftover && fileSources.Source == SourceType.OLD)
                {
                    BasePosition comparisonPosition = determineComparisonGenomePosition(
                            chromosome, position, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    position = comparisonPosition.Position;
                    chromosome = comparisonPosition.Chromosome;
                }

                int frags = var.GermlineFragments;
                int depth = (usesStart ? var.GermlineReferenceFragmentsStart : var.GermlineReferenceFragmentsEnd) + frags;
                int qual = (int)round(var.QualScore);

                BreakendData breakendData = new BreakendData(
                        breakend, var.VcfId, var.Type, chromosome, position,
                        usesStart ? var.OrientStart : var.OrientEnd, homologyOffsets, depth, frags, qual);

                geneBreakends.add(breakendData);
            }

            for(Map.Entry<String,List<BreakendData>> entry : geneBreakendMap.entrySet())
            {
                String geneName = entry.getKey();
                List<BreakendData> geneBreakends = entry.getValue();

                DisruptionData disruptionData = new DisruptionData(GERMLINE_SV, geneName, geneBreakends, mFields);
                items.add(disruptionData);
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx germline SV data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }
}
