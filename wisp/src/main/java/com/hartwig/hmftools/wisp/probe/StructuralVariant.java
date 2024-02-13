package com.hartwig.hmftools.wisp.probe;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.linx.DriverEventType.DEL;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.probe.CategoryType.AMP;
import static com.hartwig.hmftools.wisp.probe.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.FUSION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.DEFAULT_SV_BREAKENDS_PER_GENE;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.MAX_INSERT_BASES;
import static com.hartwig.hmftools.wisp.probe.ProbeConstants.MAX_POLY_A_T_BASES;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxCommonTypes;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

public class StructuralVariant extends Variant
{
    private final StructuralVariantData mVariant;
    private final List<LinxBreakend> mBreakends;
    private final List<LinxFusion> mFusions;
    private CategoryType mCategoryType;

    private List<String> mRefSequences;

    public StructuralVariant(
            final StructuralVariantData variant, final List<LinxBreakend> breakends, final List<LinxFusion> fusions)
    {
        mVariant = variant;
        mBreakends = breakends;
        mFusions = fusions;
        mRefSequences = Lists.newArrayListWithExpectedSize(2);
        mCategoryType = OTHER_SV;
        setCategoryType();
    }

    public StructuralVariantData variantData() { return mVariant; }

    public void markAmpDelDriver(boolean isAmp)
    {
        mCategoryType = isAmp ? CategoryType.AMP : CategoryType.DEL;
    }

    @Override
    public CategoryType categoryType() { return mCategoryType; }

    private void setCategoryType()
    {
        if(!mFusions.isEmpty())
        {
            mCategoryType = FUSION;
        }
        else if(mCategoryType == AMP || mCategoryType == CategoryType.DEL)
        {
            return;
        }
        else if(mBreakends.stream().anyMatch(x -> x.reportedDisruption()))
        {
            mCategoryType = DISRUPTION;
        }
        else
        {
            mCategoryType = OTHER_SV;
        }
    }

    @Override
    public String description()
    {
        if(mVariant.type() == SGL)
        {
            return format("%s %s:%d:%d",
                    mVariant.type(), mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation());
        }
        else
        {
            return format("%s %s:%d:%d - %s:%d:%d",
                    mVariant.type(), mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation());
        }
    }

    @Override
    public String gene()
    {
        if(!mFusions.isEmpty())
            return mFusions.get(0).name();

        LinxBreakend breakend = mBreakends.stream().filter(x -> x.reportedDisruption()).findFirst().orElse(null);

        if(breakend != null)
            return breakend.gene();

        return !mBreakends.isEmpty() ? mBreakends.get(0).gene() : "";
    }

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
        return mCategoryType == FUSION || mCategoryType == AMP || mCategoryType == CategoryType.DEL || mCategoryType == DISRUPTION;
    }

    @Override
    public String otherData()
    {
        return format("GcRefMin=%.2f GcRefMax=%.2f",
                mRefSequences.stream().mapToDouble(x -> calcGcPercent(x)).min().orElse(0),
                mRefSequences.stream().mapToDouble(x -> calcGcPercent(x)).max().orElse(0));
    }

    protected static List<String> generateSvReferenceSequences(
            final RefGenomeInterface refGenome, final ProbeConfig config,
            final String chrStart, final int positionStart, final String chrEnd, final int positionEnd)
    {
        List<String> refSequences = Lists.newArrayList();

        int probeLength = config.ProbeLength;
        int halfProbeLength = probeLength / 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            String chromosome = se == SE_START ? chrStart : chrEnd;
            int position = se == SE_START ? positionStart : positionEnd;

            if(position <= 0) // ie SGLs
                continue;

            int probeStart = position - halfProbeLength;

            refSequences.add(refGenome.getBaseString(chromosome, probeStart, probeStart + probeLength - 1));
        }

        return refSequences;
    }

    private static String generateSglSequence(
            final RefGenomeInterface refGenome, final ProbeConfig config,
            final String chromosome, final int position, final byte orientation, final String insertSequence)
    {
        int probeLength = config.ProbeLength;
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        // +1 take as-is
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        String sequence;

        if(orientation == POS_ORIENT)
        {
            int probeStart = position - refBaseLength + 1;
            String refBases = refGenome.getBaseString(chromosome, probeStart, position);
            sequence = refBases + insertSequence.substring(0, insSeqLength);
        }
        else
        {
            String refBases = refGenome.getBaseString(chromosome, position, position + refBaseLength - 1);
            sequence = insertSequence.substring(insertSequence.length() - insSeqLength) + refBases;
        }

        if(sequence.length() != probeLength)
        {
            CT_LOGGER.error("variant({}:{}) invalid sequenceLength({}): {}", chromosome, position, sequence.length(), sequence);
        }

        return sequence;
    }

    protected static String generateSvSequence(
            final RefGenomeInterface refGenome, final ProbeConfig config,
            final String chrStart, final int positionStart, final byte orientStart,
            final String chrEnd, final int positionEnd, final byte orientEnd, final String insertSequence)
    {
        int probeLength = config.ProbeLength;
        int halfProbeLength = probeLength / 2;
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
                basesEnd = Nucleotides.reverseComplementBases(basesEnd);
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
                basesStart = Nucleotides.reverseComplementBases(basesStart);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                basesEnd = refGenome.getBaseString(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
        }

        String sequence = basesStart + insertSequence + basesEnd;

        if(sequence.length() != probeLength)
        {
            CT_LOGGER.error("variant({}:{} - {}:{}) invalid sequenceLength({}): {}",
                    chrStart, positionStart, chrEnd, positionEnd, sequence.length(), sequence);
        }

        return sequence;
    }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome, final ProbeConfig config)
    {
        if(mCategoryType != OTHER_SV)
        {
            mRefSequences.addAll(generateSvReferenceSequences(
                    refGenome, config, mVariant.startChromosome(), mVariant.startPosition(), mVariant.endChromosome(), mVariant.endPosition()));
        }

        String sequence;

        if(mVariant.type() == SGL)
        {
            sequence = generateSglSequence(
                    refGenome, config, mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(), mVariant.insertSequence());
        }
        else
        {
            sequence = generateSvSequence(
                    refGenome, config, mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation(), mVariant.insertSequence());
        }

        setSequence(sequence);
    }

    @Override
    boolean checkFilters() { return mCategoryType != FUSION && mCategoryType != AMP && mCategoryType != CategoryType.DEL; }

    @Override
    public boolean passNonReportableFilters(final ProbeConfig config, boolean useLowerLimits)
    {
        if(reported() && mCategoryType != DISRUPTION)
            return true;

        if(!passesGcRatioLimit(gc(), config, useLowerLimits))
            return false;

        for(String refSequence : mRefSequences)
        {
            double gcRatio = calcGcPercent(refSequence);

            if(!passesGcRatioLimit(gcRatio, config, useLowerLimits))
                return false;

            if(exceedsPolyAtThreshold(refSequence))
                return false;
        }

        if(vaf() < config.VafMin)
            return false;

        if(!passesFragmentCountLimit(tumorFragments(), config, useLowerLimits))
            return false;

        return true;
    }

    private boolean exceedsPolyAtThreshold(final String sequence)
    {
        int aCount = 0;
        int tCount = 0;
        for(int i = 0; i < sequence.length(); ++i)
        {
            if(sequence.charAt(i) == 'A')
            {
                ++aCount;
            }
            if(sequence.charAt(i) == 'T')
            {
                ++tCount;
            }
            else
            {
                aCount = 0;
                tCount = 0;
            }
        }

        return aCount > MAX_POLY_A_T_BASES || tCount > MAX_POLY_A_T_BASES;
    }

    @Override
    public boolean checkAndRegisterLocation(final ProximateLocations registeredLocations)
    {
        if(registeredLocations.isNearRegisteredLocation(mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation())
        || registeredLocations.isNearRegisteredLocation(mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation()))
        {
            return false;
        }

        registeredLocations.addRegisteredLocation(mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation());
        registeredLocations.addRegisteredLocation(mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation());
        return true;
    }

    @Override
    public boolean checkAndRegisterGeneLocation(final Map<String,Integer> geneDisruptions)
    {
        for(LinxBreakend breakend : mBreakends)
        {
            Integer breakendCount = geneDisruptions.get(breakend.gene());

            if(breakendCount != null)
            {
                if(breakendCount >= DEFAULT_SV_BREAKENDS_PER_GENE)
                    return false;

                geneDisruptions.put(breakend.gene(), breakendCount + 1);
            }
            else
            {
                geneDisruptions.put(breakend.gene(), 1);
            }
        }

        return true;
    }

    public String toString()
    {
        return format("variant(%s) category(%s) fusion(%d) breakends(%d)", description(), categoryType(), mFusions.size(), mBreakends.size());
    }

    public static List<Variant> loadStructuralVariants(final String sampleId, final ProbeConfig config) throws Exception
    {
        List<Variant> variants = Lists.newArrayList();

        if(config.LinxDir == null)
            return variants;

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info
        String purpleDir = ProbeConfig.getSampleFilePath(sampleId, config.PurpleDir);

        String vcfFile = PurpleCommon.purpleSomaticSvFile(purpleDir, sampleId);

        List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(
                StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter()));

        String linxDir = ProbeConfig.getSampleFilePath(sampleId, config.LinxDir);

        List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(linxDir, sampleId));
        List<LinxSvAnnotation> annotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(linxDir, sampleId));
        List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(linxDir, sampleId));

        List<LinxDriver> drivers = LinxDriver.read(LinxDriver.generateFilename(linxDir, sampleId))
                .stream().filter(x -> x.eventType() == DEL || x.eventType() == GAIN)
                .collect(Collectors.toList());

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        if(drivers.stream().anyMatch(x -> x.eventType() == DEL))
        {
            String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, sampleId);
            GeneCopyNumberFile.read(geneCopyNumberFile).stream()
                    .filter(x -> drivers.stream().anyMatch(y -> y.gene().equals(x.geneName())))
                    .forEach(x -> geneCopyNumbers.add(x));
        }

        List<LinxCluster> clusters = LinxCluster.read(LinxCluster.generateFilename(linxDir, sampleId));

        Map<Integer,List<StructuralVariant>> clusterSVs = Maps.newHashMap();
        drivers.forEach(x -> clusterSVs.put(x.clusterId(), Lists.newArrayList()));

        for(EnrichedStructuralVariant variant : enrichedVariants)
        {
            if(variant.type() == StructuralVariantType.INF)
                continue;

            if(!variant.filter().equals(PASS))
                continue;

            if(variant.insertSequence().length() >= MAX_INSERT_BASES && variant.type() != SGL)
                continue;

            LinxSvAnnotation annotation = annotations.stream().filter(x -> x.vcfId().equals(variant.id())).findFirst().orElse(null);

            if(annotation == null)
            {
                CT_LOGGER.error("sample({}) vcfId({}) Linx annotation not found", sampleId, variant.id());
                // return Lists.newArrayList();
                continue;
            }

            List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == annotation.svId()).collect(Collectors.toList());

            List<LinxFusion> svFusions = fusions.stream()
                    .filter(x -> x.reported())
                    .filter(x -> x.chainLinks() == 0)
                    .filter(x -> svBreakends.stream().anyMatch(y -> y.id() == x.fivePrimeBreakendId()))
                    .collect(Collectors.toList());

            // only use SGLs if in a reportable fusion
            if(variant.type() == SGL && svFusions.isEmpty())
                 continue;

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

            if(cluster == null || cluster.category().equals(LinxCommonTypes.SUPER_TYPE_ARTIFACT))
                continue;

            StructuralVariantData variantData = convertSvData(variant, annotation.svId());

            StructuralVariant sv = new StructuralVariant(variantData, svBreakends, svFusions);
            variants.add(sv);

            if(clusterSVs.containsKey(cluster.clusterId()))
            {
                clusterSVs.get(cluster.clusterId()).add(sv);
            }
        }

        CT_LOGGER.info("loaded {} structural variants from vcf({})", variants.size(), vcfFile);

        // find SVs related to DEL and AMP events
        for(LinxDriver driver : drivers)
        {
            List<StructuralVariant> svList = clusterSVs.get(driver.clusterId());
            StructuralVariant driverSv = null;

            if(driver.eventType() == GAIN)
            {
                double maxJcn = 0;
                for(StructuralVariant sv : svList)
                {
                    double svJcn = max(sv.variantData().adjustedStartCopyNumberChange(), sv.variantData().adjustedEndCopyNumberChange());
                    if(svJcn > maxJcn)
                    {
                        maxJcn = svJcn;
                        driverSv = sv;
                    }
                }
            }
            else
            {
                GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream().filter(x -> x.geneName().equals(driver.gene())).findFirst().orElse(null);

                if(geneCopyNumber != null)
                {
                    for(StructuralVariant sv : svList)
                    {
                        if(matchesDelRegion(sv, geneCopyNumber))
                            sv.markAmpDelDriver(false);
                    }
                }
            }

            if(driverSv != null)
                driverSv.markAmpDelDriver(true);
        }

        return variants;
    }

    public static boolean matchesDelRegion(final StructuralVariant sv, final GeneCopyNumber geneCopyNumber)
    {
        if(sv.variantData().startOrientation() == POS_ORIENT && abs(sv.variantData().startPosition() - geneCopyNumber.minRegionStart()) <= 1)
            return true;

        if(sv.variantData().endOrientation() == POS_ORIENT && abs(sv.variantData().endPosition() - geneCopyNumber.minRegionStart()) <= 1)
            return true;

        if(sv.variantData().startOrientation() == NEG_ORIENT && abs(sv.variantData().startPosition() - geneCopyNumber.minRegionEnd()) <= 1)
            return true;

        if(sv.variantData().endOrientation() == NEG_ORIENT && abs(sv.variantData().endPosition() - geneCopyNumber.minRegionEnd()) <= 1)
            return true;

        return false;
    }
}
