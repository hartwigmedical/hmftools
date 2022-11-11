package com.hartwig.hmftools.ctdna;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.ctdna.CategoryType.FUSION;
import static com.hartwig.hmftools.ctdna.CategoryType.GERMLINE_SV;
import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.ctdna.PvConfig.MAX_INSERT_BASES;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.StructuralVariant.generateSvReferenceSequences;
import static com.hartwig.hmftools.ctdna.StructuralVariant.generateSvSequence;
import static com.hartwig.hmftools.ctdna.VariantSelection.addRegisteredLocation;
import static com.hartwig.hmftools.ctdna.VariantSelection.isNearRegisteredLocation;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxCommonTypes;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

public class GermlineSv extends Variant
{
    private final LinxGermlineSv mVariant;

    private List<String> mRefSequences;

    public GermlineSv(final LinxGermlineSv variant)
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
    public void generateSequences(final RefGenomeInterface refGenome, final PvConfig config)
    {
        mRefSequences.addAll(generateSvReferenceSequences(
                refGenome, config, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.ChromosomeEnd, mVariant.PositionEnd));

        String sequence = generateSvSequence(
                refGenome, config, mVariant.ChromosomeStart, mVariant.PositionStart, mVariant.OrientStart,
                mVariant.ChromosomeEnd, mVariant.PositionEnd, mVariant.OrientEnd, mVariant.InsertSequence);

        setSequence(sequence);
    }

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

    public static List<Variant> loadGermlineStructuralVariants(final String sampleId, final PvConfig config) throws Exception
    {
        List<Variant> variants = Lists.newArrayList();

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info
        String filename = LinxGermlineSv.generateFilename(config.LinxDir, sampleId);
        final List<LinxGermlineSv> germlineSvs = LinxGermlineSv.read(filename);

        germlineSvs.stream().filter(x -> x.Reported).forEach(x -> variants.add(new GermlineSv(x)));

        PV_LOGGER.info("loaded {} germline SVs from vcf({})", variants.size(), filename);

        return variants;
    }
}
