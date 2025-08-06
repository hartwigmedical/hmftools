package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.max;
import static java.lang.String.format;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSvProbe;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GermlineSv extends Variant
{
    private final LinxGermlineDisruption mVariant;

    private static final Logger LOGGER = LogManager.getLogger(GermlineSv.class);

    public GermlineSv(final LinxGermlineDisruption variant)
    {
        mVariant = variant;
    }

    @Override
    public CategoryType categoryType()
    {
        return CategoryType.GERMLINE_SV;
    }

    @Override
    public String description()
    {
        return format("%s %s:%d:%d - %s:%d:%d",
                mVariant.Type, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd);
    }

    @Override
    public String gene()
    {
        return mVariant.GeneName;
    }

    @Override
    public double copyNumber()
    {
        return 0;
    }

    @Override
    public double vaf()
    {
        double refFrags = (mVariant.GermlineReferenceFragmentsStart + mVariant.GermlineReferenceFragmentsEnd) / 2.0;
        return refFrags > 0 ? mVariant.GermlineFragments / refFrags : 0;
    }

    @Override
    public int tumorFragments()
    {
        return max(mVariant.TumorReferenceFragmentsStart, mVariant.TumorReferenceFragmentsEnd);
    }

    @Override
    public boolean reported()
    {
        return true;
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
    public boolean checkFilters()
    {
        return false;
    }

    @Override
    public boolean checkAndRegisterLocation(ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart)
                || registeredLocations.isNearRegisteredLocation(mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd))
        {
            return false;
        }

        registeredLocations.addRegisteredLocation(mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart);
        registeredLocations.addRegisteredLocation(mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd);
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s)", description(), categoryType());
    }

    public static List<GermlineSv> load(final String sampleId, final String linxGermlineDir)
    {
        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info

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

        ArrayList<GermlineSv> variants = new ArrayList<>();

        List<LinxGermlineDisruption> germlineSvs;
        List<LinxBreakend> germlineBreakends;

        try
        {
            germlineSvs = LinxGermlineDisruption.read(germlineSvFile);
            germlineBreakends = LinxBreakend.read(germlineBreakendsFile).stream()
                    .filter(LinxBreakend::reportedDisruption).toList();
        }
        catch(IOException e)
        {
            String error = "Failed to load germline structural variants: " + e;
            LOGGER.error(error);
            throw new RuntimeException(error);
        }

        for(LinxGermlineDisruption germlineSv : germlineSvs)
        {
            if(germlineBreakends.stream().anyMatch(breakend -> breakend.svId() == germlineSv.SvId))
            {
                variants.add(new GermlineSv(germlineSv));
            }
        }

        LOGGER.info("loaded {} germline SVs", variants.size());

        return variants;
    }
}
