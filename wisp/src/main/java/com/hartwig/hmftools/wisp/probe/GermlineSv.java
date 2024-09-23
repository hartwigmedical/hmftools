package com.hartwig.hmftools.wisp.probe;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.probe.CategoryType.GERMLINE_SV;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;

public class GermlineSv extends Variant
{
    private final LinxGermlineDisruption mVariant;

    private List<String> mRefSequences;

    public GermlineSv(final LinxGermlineDisruption variant)
    {
        mVariant = variant;
        mRefSequences = Lists.newArrayListWithExpectedSize(2);
    }

    @Override
    public CategoryType categoryType() { return GERMLINE_SV; }

    @Override
    public String description()
    {
        return format("%s %s:%d:%d - %s:%d:%d",
                mVariant.Type, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd);
    }

    @Override
    public String gene() { return mVariant.GeneName; }

    @Override
    public List<String> refSequences() { return mRefSequences; }

    @Override
    public double copyNumber() { return 0; }

    @Override
    public double vaf()
    {
        double refFrags = (mVariant.GermlineReferenceFragmentsStart + mVariant.GermlineReferenceFragmentsEnd) / 2.0;
        return refFrags > 0 ? mVariant.GermlineFragments / refFrags : 0;
    }

    @Override
    public int tumorFragments() { return max(mVariant.TumorReferenceFragmentsStart, mVariant.TumorReferenceFragmentsEnd); }

    @Override
    public boolean hasPhaseVariants() { return false; }

    @Override
    public boolean reported() { return true; }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final ProbeConfig config)
    {
        mRefSequences.addAll(StructuralVariant.generateSvReferenceSequences(
                refGenome, config, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.ChromosomeEnd, mVariant.PositionEnd));

        String sequence = StructuralVariant.generateSvSequence(
                refGenome, config, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd, mVariant.InsertSequence);

        setSequence(sequence);
    }

    @Override
    boolean checkFilters() { return false; }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
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

    public static List<Variant> loadGermlineStructuralVariants(final String sampleId, final ProbeConfig config) throws Exception
    {
        List<Variant> variants = Lists.newArrayList();

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info
        String linxDir = ProbeConfig.getSampleFilePath(sampleId, config.LinxGermlineDir);

        if(linxDir == null)
            return variants;

        String germlineSvFile = LinxGermlineDisruption.generateFilename(linxDir, sampleId);
        String germlineBreakendsFile = LinxBreakend.generateFilename(linxDir, sampleId, true);

        if(!Files.exists(Paths.get(germlineSvFile)) || !Files.exists(Paths.get(germlineBreakendsFile)))
            return variants;

        List<LinxGermlineDisruption> germlineSvs = LinxGermlineDisruption.read(germlineSvFile);
        List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendsFile).stream()
                .filter(x -> x.reportedDisruption()).collect(Collectors.toList());

        for(LinxGermlineDisruption germlineSv : germlineSvs)
        {
            if(germlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
                variants.add(new GermlineSv(germlineSv));
        }

        CT_LOGGER.info("loaded {} germline SVs from vcf({})", variants.size(), germlineSvFile);

        return variants;
    }
}
