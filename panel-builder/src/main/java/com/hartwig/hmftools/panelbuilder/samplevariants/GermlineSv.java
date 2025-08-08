package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSvProbe;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Germline structural variant.
public class GermlineSv extends Variant
{
    private final LinxGermlineDisruption mVariant;
    private final List<LinxBreakend> mBreakends;

    private static final Logger LOGGER = LogManager.getLogger(GermlineSv.class);

    private GermlineSv(final LinxGermlineDisruption variant, final List<LinxBreakend> breakends)
    {
        mVariant = variant;
        mBreakends = breakends;
    }

    // TODO: only select if driver = true
    @Override
    public boolean isDriver()
    {
        return mBreakends.stream().anyMatch(LinxBreakend::reportedDisruption);
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        return buildSvProbe(
                mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd, mVariant.InsertSequence,
                PROBE_LENGTH, refGenome);
    }

    @Override
    public List<ProximateLocations.Location> checkedLocations()
    {
        return List.of(
                new ProximateLocations.Location(mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart),
                new ProximateLocations.Location(mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd));
    }

    @Override
    public String toString()
    {
        return format("%s %s:%d:%d - %s:%d:%d",
                mVariant.Type, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd);
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

        LOGGER.info("Loaded {} germline structural variants", variants.size());
        variants.forEach(variant -> LOGGER.trace("GermlineSv: {}", variant));

        return variants;
    }
}
