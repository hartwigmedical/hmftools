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
import static com.hartwig.hmftools.ctdna.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.ctdna.PvConfig.MAX_INSERT_BASES;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
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
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

public class StructuralVariant implements Variant
{
    private final StructuralVariantData mVariant;
    private final List<LinxBreakend> mBreakends;
    private final List<LinxFusion> mFusions;
    private final LinxCluster mCluster;

    private String mSequence;
    private List<String> mRefSequences;

    public StructuralVariant(
            final StructuralVariantData variant, final List<LinxBreakend> breakends, final List<LinxFusion> fusions, final LinxCluster cluster)
    {
        mVariant = variant;
        mBreakends = breakends;
        mFusions = fusions;
        mCluster = cluster;
        mSequence = "";
        mRefSequences = Lists.newArrayListWithExpectedSize(2);
    }

    @Override
    public CategoryType categoryType()
    {
        if(mFusions.stream().anyMatch(x -> x.reported()))
            return FUSION;

        return OTHER_SV;
    }

    @Override
    public String description()
    {
        return format("%s %s:%d%d - %s:%d:%d",
                mVariant.type(), mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation());
    }

    @Override
    public String gene()
    {
        LinxFusion fusion = mFusions.stream().filter(x -> x.reported()).findFirst().orElse(null);

        if(fusion != null)
            return fusion.name();

        LinxBreakend breakend = mBreakends.stream().filter(x -> x.reportedDisruption()).findFirst().orElse(null);

        if(breakend != null)
            return breakend.gene();

        return !mBreakends.isEmpty() ? mBreakends.get(0).gene() : "";
    }

    @Override
    public String sequence() { return mSequence; }

    @Override
    public List<String> refSequences() { return mRefSequences; }

    @Override
    public double copyNumber() { return max(mVariant.adjustedStartCopyNumberChange(), mVariant.adjustedEndCopyNumberChange()); }

    @Override
    public double vaf() { return mVariant.adjustedStartAF(); }

    @Override
    public int tumorFragments() { return max(mVariant.startTumorVariantFragmentCount(), mVariant.endTumorVariantFragmentCount()); }

    @Override
    public boolean hasPhaseVariants() { return false; }

    @Override
    public boolean reported()
    {
        if(mFusions.stream().anyMatch(x -> x.reported()))
            return true;

        if(mBreakends.stream().anyMatch(x -> x.reportedDisruption()))
            return true;

        return false;
    }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final PvConfig config)
    {
        int positionStart = mVariant.startPosition();
        byte orientStart = mVariant.startOrientation();
        String chrStart = mVariant.startChromosome();

        int positionEnd = mVariant.endPosition();
        byte orientEnd = mVariant.endOrientation();
        String chrEnd = mVariant.endChromosome();

        int probeLength = config.ProbeLength;
        int halfProbeLength = probeLength / 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            String chromosome = se == SE_START ? chrStart : chrEnd;
            int position = se == SE_START ? positionStart : positionEnd;
            int probeStart = position - halfProbeLength;

            mRefSequences.add(refGenome.getBaseString(chromosome, probeStart, probeStart + probeLength - 1));
        }

        String insertSequence = mVariant.insertSequence();
        int insSeqLength = insertSequence.length();
        int halfInsSeqLength = insSeqLength / 2;
        int halfNonInsSeqLength = halfProbeLength - halfInsSeqLength;

        // +1 to -1 - do nothing
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        String basesStart;
        String basesEnd;

        if(orientStart == POS_ORIENT)
        {
            int probeStart = positionStart - halfNonInsSeqLength + 1;
            basesStart = refGenome.getBaseString(chrStart, probeStart, positionStart);

            int endBaseLength = probeLength - basesStart.length() - insSeqLength;

            if(orientEnd == NEG_ORIENT)
            {
                basesEnd = refGenome.getBaseString(chrEnd, positionEnd, positionEnd + endBaseLength - 1);
            }
            else
            {
                basesEnd = refGenome.getBaseString(chrEnd, positionEnd - endBaseLength + 1, positionEnd);
                basesEnd = Nucleotides.reverseStrandBases(basesEnd);
            }
        }
        else
        {
            if(orientEnd == POS_ORIENT)
            {
                // swap ends and treat as +1/-1
                int probeStart = positionEnd - halfNonInsSeqLength + 1;
                basesStart = refGenome.getBaseString(chrEnd, probeStart, positionEnd);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                basesEnd = refGenome.getBaseString(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                basesStart = refGenome.getBaseString(chrEnd, positionEnd, positionEnd + halfNonInsSeqLength - 1);
                basesStart = Nucleotides.reverseStrandBases(basesStart);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                basesEnd = refGenome.getBaseString(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
        }

        mSequence = basesStart + insertSequence + basesEnd;

        if(mSequence.length() != probeLength)
        {
            PV_LOGGER.error("variant({}) invalid sequenceLength({}): {}", description(), mSequence.length(), mSequence);
            return;
        }
    }

    @Override
    public double gc() { return VariantUtils.calcGcPercent(mSequence); }

    @Override
    public boolean checkAndRegisterLocation(final Map<String,List<Integer>> registeredLocations)
    {
        if(isNearRegisteredLocation(registeredLocations, mVariant.startChromosome(), mVariant.startPosition())
        || isNearRegisteredLocation(registeredLocations, mVariant.endChromosome(), mVariant.endPosition()))
        {
            return false;
        }

        addRegisteredLocation(registeredLocations, mVariant.startChromosome(), mVariant.startPosition());
        addRegisteredLocation(registeredLocations, mVariant.endChromosome(), mVariant.endPosition());
        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s) fusion(%d) breakends(%d)", description(), categoryType(), mFusions.size(), mBreakends.size());
    }

    public static List<Variant> loadStructuralVariants(final String sampleId, final PvConfig config)
    {
        List<Variant> variants = Lists.newArrayList();

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info
        try
        {
            String purpleDir = PvConfig.getSampleFilePath(sampleId, config.PurpleDir);
            String linxDir = PvConfig.getSampleFilePath(sampleId, config.LinxDir);

            String vcfFile = PurpleCommon.purpleSvFile(purpleDir, sampleId);

            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(
                    StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter()));

            List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(linxDir, sampleId));
            List<LinxSvAnnotation> annotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(linxDir, sampleId));
            List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(linxDir, sampleId));

            List<LinxCluster> clusters = LinxCluster.read(LinxCluster.generateFilename(linxDir, sampleId));

            for(EnrichedStructuralVariant variant : enrichedVariants)
            {
                if(variant.type() == StructuralVariantType.INF || variant.type() == StructuralVariantType.SGL)
                    continue;

                if(!variant.filter().equals(PASS))
                    continue;

                if(variant.insertSequence().length() >= MAX_INSERT_BASES)
                    continue;

                LinxSvAnnotation annotation = annotations.stream().filter(x -> x.vcfId().equals(variant.id())).findFirst().orElse(null);

                if(annotation == null)
                {
                    PV_LOGGER.error("sample({}) vcfId({}) Linx annotation not found", sampleId, variant.id());
                    return Lists.newArrayList();
                }

                List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == annotation.svId()).collect(Collectors.toList());

                List<LinxFusion> svFusions = fusions.stream()
                        .filter(x -> x.chainLinks() == 0)
                        .filter(x -> svBreakends.stream().anyMatch(y -> y.id() == x.fivePrimeBreakendId()))
                        .collect(Collectors.toList());

                LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

                if(cluster == null || cluster.category().equals(LinxCommonTypes.SUPER_TYPE_ARTIFACT))
                    continue;

                StructuralVariantData variantData = convertSvData(variant, annotation.svId());

                variants.add(new StructuralVariant(variantData, svBreakends, svFusions, cluster));
            }

            PV_LOGGER.info("loaded {} structural variants from vcf({})", variants.size(), vcfFile);
        }
        catch(Exception e)
        {
            PV_LOGGER.error("sample({}) failed to load Purple or Linx SV files: {}", sampleId, e.toString());
        }

        return variants;
    }
}
