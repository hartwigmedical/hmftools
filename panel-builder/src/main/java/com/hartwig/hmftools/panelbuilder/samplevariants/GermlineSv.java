package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSvProbe;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.panelbuilder.SequenceDefinition;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Germline structural variant.
public class GermlineSv implements Variant
{
    private final LinxGermlineDisruption mVariant;
    private final List<LinxBreakend> mBreakends;

    private static final Logger LOGGER = LogManager.getLogger(GermlineSv.class);

    private GermlineSv(final LinxGermlineDisruption variant, final List<LinxBreakend> breakends)
    {
        mVariant = variant;
        mBreakends = breakends;
    }

    private String gene()
    {
        return mVariant.GeneName;
    }

    @Override
    public boolean isDriver()
    {
        return mBreakends.stream().anyMatch(LinxBreakend::reportedDisruption);
    }

    @Override
    public SequenceDefinition generateProbe()
    {
        return buildSvProbe(
                mVariant.ChromosomeStart, mVariant.PositionStart, Orientation.fromByte(mVariant.OrientStart),
                mVariant.ChromosomeEnd, mVariant.PositionEnd, Orientation.fromByte(mVariant.OrientEnd), mVariant.InsertSequence,
                PROBE_LENGTH);
    }

    @Override
    public TargetMetadata.Type targetType()
    {
        if(isDriver())
        {
            return TargetMetadata.Type.SAMPLE_GERMLINE_SV_DRIVER;
        }
        else
        {
            // Shouldn't happen because nondrivers are filtered out.
            throw new IllegalStateException("Unhandled germline SV type");
        }
    }

    @Override
    public String toString()
    {
        return format("%s:%d:%d - %s:%d:%d %s %s",
                mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd,
                mVariant.Type, Optional.ofNullable(gene()).orElse(""));
    }

    public static List<GermlineSv> load(final String sampleId, final String linxGermlineDir)
    {
        if(linxGermlineDir == null)
        {
            return emptyList();
        }

        String germlineSvFile = LinxGermlineDisruption.generateFilename(linxGermlineDir, sampleId);
        String germlineBreakendsFile = LinxBreakend.generateFilename(linxGermlineDir, sampleId, true);

        if(!Files.exists(Paths.get(germlineSvFile)) || !Files.exists(Paths.get(germlineBreakendsFile)))
        {
            return emptyList();
        }

        List<LinxGermlineDisruption> germlineSvs;
        List<LinxBreakend> germlineBreakends;
        try
        {
            germlineSvs = LinxGermlineDisruption.read(germlineSvFile);
            germlineBreakends = LinxBreakend.read(germlineBreakendsFile);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to load germline structural variants: " + e);
        }

        List<GermlineSv> variants = germlineSvs.stream()
                .map(germlineSv ->
                {
                    List<LinxBreakend> svBreakends = germlineBreakends.stream()
                            .filter(breakend -> breakend.svId() == germlineSv.SvId)
                            .toList();
                    return new GermlineSv(germlineSv, svBreakends);
                })
                .toList();

        LOGGER.debug("Loaded {} germline structural variants", variants.size());
        variants.forEach(variant -> LOGGER.trace("GermlineSv: {}", variant));

        return variants;
    }
}
