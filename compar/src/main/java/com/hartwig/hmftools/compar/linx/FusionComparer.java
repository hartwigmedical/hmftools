package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.linx.DisruptionComparer.buildBreakendData;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class FusionComparer extends ItemComparer
{
    protected enum Fields
    {
        Reported,
        ReportedType,
        Phased,
        Likelihood,
        FusedExonUp,
        FusedExonDown,
        ChainLinks,
        ChainTerminated,
        DomainsKept,
        DomainsLost,
        BreakendUp,
        BreakendDown;
    }

    private DisruptionComparer mDisruptionComparer;

    public FusionComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Reported.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Reported.toString()), null));

        mFields.add(new FieldInfo(
                Fields.ReportedType.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.ReportedType.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Phased.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Phased.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Likelihood.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Likelihood.toString()), null));

        mFields.add(new FieldInfo(
                Fields.FusedExonUp.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FusedExonUp.toString(), null, null), null));

        mFields.add(new FieldInfo(
                Fields.FusedExonDown.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.FusedExonDown.toString(), null, null), null));

        mFields.add(new FieldInfo(
                Fields.ChainLinks.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.ChainLinks.toString(), null, null), null));

        mFields.add(new FieldInfo(
                Fields.ChainTerminated.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.ChainTerminated.toString()), null));

        mFields.add(new FieldInfo(
                Fields.DomainsKept.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DomainsKept.toString()), null));

        mFields.add(new FieldInfo(
                Fields.DomainsLost.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.DomainsLost.toString()), null));

        // TODO: display only
        mFields.add(new FieldInfo(
                Fields.BreakendUp.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.BreakendUp.toString()), null));

        mFields.add(new FieldInfo(
                Fields.BreakendDown.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.BreakendDown.toString()), null));

        mDisruptionComparer = null;
    }

    public void setDisruptionComparer(final DisruptionComparer disruptionComparer) { mDisruptionComparer = disruptionComparer; }

    @Override
    public CategoryType category() { return FUSION; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());

        // excluded unless matching breakends can be loaded: FLD_TRANSCRIPT_UP, FLD_TRANSCRIPT_DOWN, FLD_JUNCTION_COPY_NUMBER
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        try
        {
            List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(fileSources.Linx, sampleId));
            return processFusions(fusions, fileSources.Source);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx fusion data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> processFusions(final List<LinxFusion> fusions, final SourceType sourceType)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        List<LinxBreakend> breakends = mDisruptionComparer != null ?
                mDisruptionComparer.breakends().get(sourceType) : Collections.emptyList();


        List<StructuralVariantData> svDataList = mDisruptionComparer != null ?
                mDisruptionComparer.svDataList().get(sourceType) : Collections.emptyList();

        for(LinxFusion fusion : fusions)
        {
            if(fusion.reportedType().equals(KnownFusionType.NONE.toString()))
                continue;

            String[] geneNames = LinxFusion.geneNames(fusion.name());

            BreakendData breakendStart = buildBreakend(fusion.fivePrimeBreakendId(), geneNames[0], breakends, svDataList);
            BreakendData breakendEnd = buildBreakend(fusion.threePrimeBreakendId(), geneNames[1], breakends, svDataList);

            comparableItems.add(new FusionData(fusion, fusion.name(), breakendStart, breakendEnd, mFields));
        }

        return comparableItems;
    }

    private BreakendData buildBreakend(
            final int breakendId, final String geneName, final List<LinxBreakend> breakends, final List<StructuralVariantData> svDataList)
    {
        LinxBreakend breakend = breakends.stream()
                .filter(x -> x.id() == breakendId && x.gene().equals(geneName)).findFirst().orElse(null);

        if(breakend == null)
        {
            breakend = breakends.stream().filter(x -> x.id() == breakendId).findFirst().orElse(null);

            if(breakend == null)
                return null;
        }

        int svId = breakend.svId();

        StructuralVariantData svData = svDataList.stream().filter(x -> x.id() == svId).findFirst().orElse(null);

        if(svData == null)
            return null;

        return buildBreakendData(breakend, svData, null);
    }
}
