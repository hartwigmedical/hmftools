package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.wisp.CategoryType.GERMLINE_SV;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeGenerator.generateSvProbe;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.wisp.CategoryType;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.ProbeFactory;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;

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
        return GERMLINE_SV;
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
    public void generateProbe(final RefGenomeInterface refGenome, final ProbeFactory probeFactory, final PanelCoverage coverage)
    {
        ProbeGenerationResult result = generateSvProbe(
                mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd, mVariant.InsertSequence,
                targetMetadata(), refGenome, probeFactory, coverage);
        setProbeGenResult(result);
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

    public static List<Variant> loadGermlineStructuralVariants(final String sampleId, final String linxGermlineDir)
    {
        List<Variant> variants = Lists.newArrayList();

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info

        if(linxGermlineDir == null)
        {
            return variants;
        }

        String germlineSvFile = LinxGermlineDisruption.generateFilename(linxGermlineDir, sampleId);
        String germlineBreakendsFile = LinxBreakend.generateFilename(linxGermlineDir, sampleId, true);

        if(!Files.exists(Paths.get(germlineSvFile)) || !Files.exists(Paths.get(germlineBreakendsFile)))
        {
            return variants;
        }

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
            if(germlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
            {
                variants.add(new GermlineSv(germlineSv));
            }
        }

        LOGGER.info("loaded {} germline SVs from vcf({})", variants.size(), germlineSvFile);

        return variants;
    }
}
