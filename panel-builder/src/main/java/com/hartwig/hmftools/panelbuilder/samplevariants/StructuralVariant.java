package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.linx.DriverEventType.DEL;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.wisp.CategoryType.AMP;
import static com.hartwig.hmftools.common.wisp.CategoryType.DISRUPTION;
import static com.hartwig.hmftools.common.wisp.CategoryType.FUSION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.MAX_INSERT_BASES;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.MAX_POLY_A_T_BASES;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.SV_BREAKENDS_PER_GENE;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.VAF_MIN;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
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
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class StructuralVariant extends Variant
{
    private final StructuralVariantData mVariant;
    private final List<LinxBreakend> mBreakends;
    private final List<LinxFusion> mFusions;
    private CategoryType mCategoryType;

    private final List<String> mRefSequences;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariant.class);

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

    private StructuralVariantData variantData()
    {
        return mVariant;
    }

    private void markAmpDelDriver(boolean isAmp)
    {
        mCategoryType = isAmp ? CategoryType.AMP : CategoryType.DEL;
    }

    @Override
    public CategoryType categoryType()
    {
        return mCategoryType;
    }

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
        else if(mBreakends.stream().anyMatch(LinxBreakend::reportedDisruption))
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
        {
            return mFusions.get(0).name();
        }

        LinxBreakend breakend = mBreakends.stream().filter(LinxBreakend::reportedDisruption).findFirst().orElse(null);

        if(breakend != null)
        {
            return breakend.gene();
        }

        return !mBreakends.isEmpty() ? mBreakends.get(0).gene() : "";
    }

    @Override
    public List<String> refSequences()
    {
        return mRefSequences;
    }

    @Override
    public double copyNumber()
    {
        return max(mVariant.adjustedStartCopyNumberChange(), mVariant.adjustedEndCopyNumberChange());
    }

    @Override
    public double vaf()
    {
        return mVariant.adjustedStartAF();
    }

    @Override
    public int tumorFragments()
    {
        return max(mVariant.startTumorVariantFragmentCount(), mVariant.endTumorVariantFragmentCount());
    }

    @Override
    public boolean reported()
    {
        return mCategoryType == FUSION || mCategoryType == AMP || mCategoryType == CategoryType.DEL || mCategoryType == DISRUPTION;
    }

    @Override
    public String otherData()
    {
        return format("GcRefMin=%.2f GcRefMax=%.2f",
                mRefSequences.stream().mapToDouble(GcCalcs::calcGcPercent).min().orElse(0),
                mRefSequences.stream().mapToDouble(GcCalcs::calcGcPercent).max().orElse(0));
    }

    protected static List<String> generateSvReferenceSequences(
            final RefGenomeInterface refGenome,
            final String chrStart, final int positionStart, final String chrEnd, final int positionEnd)
    {
        List<String> refSequences = Lists.newArrayList();

        int probeLength = PROBE_LENGTH;
        int halfProbeLength = probeLength / 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            String chromosome = se == SE_START ? chrStart : chrEnd;
            int position = se == SE_START ? positionStart : positionEnd;

            if(position <= 0) // ie SGLs
            {
                continue;
            }

            int probeStart = position - halfProbeLength;

            refSequences.add(refGenome.getBaseString(chromosome, probeStart, probeStart + probeLength - 1));
        }

        return refSequences;
    }

    private static String generateSglSequence(
            final RefGenomeInterface refGenome,
            final String chromosome, final int position, final byte orientation, final String insertSequence)
    {
        int probeLength = PROBE_LENGTH;
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        // +1 take as-is
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        String sequence;

        if(orientation == ORIENT_FWD)
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
            LOGGER.error("variant({}:{}) invalid sequenceLength({}): {}", chromosome, position, sequence.length(), sequence);
        }

        return sequence;
    }

    protected static String generateSvSequence(
            final RefGenomeInterface refGenome,
            final String chrStart, final int positionStart, final byte orientStart,
            final String chrEnd, final int positionEnd, final byte orientEnd, final String insertSequence)
    {
        int probeLength = PROBE_LENGTH;
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

        if(orientStart == ORIENT_FWD)
        {
            int probeStart = positionStart - halfNonInsSeqLength + 1;
            basesStart = refGenome.getBaseString(chrStart, probeStart, positionStart);

            int endBaseLength = probeLength - basesStart.length() - insSeqLength;

            if(orientEnd == ORIENT_REV)
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
            if(orientEnd == ORIENT_FWD)
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
            LOGGER.error("variant({}:{} - {}:{}) invalid sequenceLength({}): {}",
                    chrStart, positionStart, chrEnd, positionEnd, sequence.length(), sequence);
        }

        return sequence;
    }

    @Override
    public void generateSequences(final RefGenomeInterface refGenome)
    {
        if(mCategoryType != OTHER_SV)
        {
            mRefSequences.addAll(generateSvReferenceSequences(
                    refGenome, mVariant.startChromosome(), mVariant.startPosition(), mVariant.endChromosome(), mVariant.endPosition()));
        }

        String sequence;

        if(mVariant.type() == SGL)
        {
            sequence = generateSglSequence(
                    refGenome, mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(), mVariant.insertSequence());
        }
        else
        {
            sequence = generateSvSequence(
                    refGenome, mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation(), mVariant.insertSequence());
        }

        setSequence(sequence);
    }

    @Override
    public boolean checkFilters()
    {
        return mCategoryType != FUSION && mCategoryType != AMP && mCategoryType != CategoryType.DEL;
    }

    @Override
    public boolean passNonReportableFilters(boolean useLowerLimits)
    {
        if(reported() && mCategoryType != DISRUPTION)
        {
            return true;
        }

        if(!passesGcRatioLimit(gc(), useLowerLimits))
        {
            return false;
        }

        for(String refSequence : mRefSequences)
        {
            double gcRatio = calcGcPercent(refSequence);

            if(!passesGcRatioLimit(gcRatio, useLowerLimits))
            {
                return false;
            }

            if(exceedsPolyAtThreshold(refSequence))
            {
                return false;
            }
        }

        if(vaf() < VAF_MIN)
        {
            return false;
        }

        if(!passesFragmentCountLimit(tumorFragments(), useLowerLimits))
        {
            return false;
        }

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
    public boolean checkAndRegisterGeneLocation(final Map<String, Integer> geneDisruptions)
    {
        for(LinxBreakend breakend : mBreakends)
        {
            Integer breakendCount = geneDisruptions.get(breakend.gene());

            if(breakendCount != null)
            {
                if(breakendCount >= SV_BREAKENDS_PER_GENE)
                {
                    return false;
                }

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

    public static List<Variant> loadStructuralVariants(final String sampleId, final String purpleDir, @Nullable final String linxDir)
    {
        List<Variant> variants = Lists.newArrayList();

        if(linxDir == null)
        {
            return variants;
        }

        // load each structural variant (ignoring INFs and SGLs), and link to any disruption/breakend and fusion, and cluster info

        String vcfFile = PurpleCommon.purpleSomaticSvFile(purpleDir, sampleId);

        List<EnrichedStructuralVariant> enrichedVariants;
        List<LinxBreakend> breakends;
        List<LinxSvAnnotation> annotations;
        List<LinxFusion> fusions;
        List<LinxDriver> drivers;
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<LinxCluster> clusters;

        try
        {
            enrichedVariants = new EnrichedStructuralVariantFactory().enrich(
                    StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter()));
            breakends = LinxBreakend.read(LinxBreakend.generateFilename(linxDir, sampleId));
            annotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(linxDir, sampleId));
            fusions = LinxFusion.read(LinxFusion.generateFilename(linxDir, sampleId));

            drivers = LinxDriver.read(LinxDriver.generateFilename(linxDir, sampleId))
                    .stream().filter(x -> x.eventType() == DEL || x.eventType() == GAIN)
                    .toList();

            if(drivers.stream().anyMatch(x -> x.eventType() == DEL))
            {
                String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, sampleId);
                GeneCopyNumberFile.read(geneCopyNumberFile).stream()
                        .filter(x -> drivers.stream().anyMatch(y -> y.gene().equals(x.geneName())))
                        .forEach(geneCopyNumbers::add);
            }

            clusters = LinxCluster.read(LinxCluster.generateFilename(linxDir, sampleId));
        }
        catch(IOException e)
        {
            String error = "Failed to load structural variants: " + e;
            LOGGER.error(error);
            throw new RuntimeException(error);
        }

        Map<Integer, List<StructuralVariant>> clusterSVs = Maps.newHashMap();
        drivers.forEach(x -> clusterSVs.put(x.clusterId(), Lists.newArrayList()));

        for(EnrichedStructuralVariant variant : enrichedVariants)
        {
            if(variant.type() == StructuralVariantType.INF)
            {
                continue;
            }

            if(!variant.filter().equals(PASS))
            {
                continue;
            }

            if(variant.insertSequence().length() >= MAX_INSERT_BASES && variant.type() != SGL)
            {
                continue;
            }

            LinxSvAnnotation annotation = annotations.stream().filter(x -> x.vcfId().equals(variant.id())).findFirst().orElse(null);

            if(annotation == null)
            {
                LOGGER.error("sample({}) vcfId({}) Linx annotation not found", sampleId, variant.id());
                // return Lists.newArrayList();
                continue;
            }

            List<LinxBreakend> svBreakends = breakends.stream().filter(x -> x.svId() == annotation.svId()).collect(Collectors.toList());

            List<LinxFusion> svFusions = fusions.stream()
                    .filter(LinxFusion::reported)
                    .filter(x -> x.chainLinks() == 0)
                    .filter(x -> svBreakends.stream().anyMatch(y -> y.id() == x.fivePrimeBreakendId()))
                    .collect(Collectors.toList());

            // only use SGLs if in a reportable fusion
            if(variant.type() == SGL && svFusions.isEmpty())
            {
                continue;
            }

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

            if(cluster == null || cluster.category().equals(LinxCommonTypes.SUPER_TYPE_ARTIFACT))
            {
                continue;
            }

            StructuralVariantData variantData = convertSvData(variant, annotation.svId());

            StructuralVariant sv = new StructuralVariant(variantData, svBreakends, svFusions);
            variants.add(sv);

            if(clusterSVs.containsKey(cluster.clusterId()))
            {
                clusterSVs.get(cluster.clusterId()).add(sv);
            }
        }

        LOGGER.info("loaded {} structural variants from vcf({})", variants.size(), vcfFile);

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
                GeneCopyNumber geneCopyNumber =
                        geneCopyNumbers.stream().filter(x -> x.geneName().equals(driver.gene())).findFirst().orElse(null);

                if(geneCopyNumber != null)
                {
                    for(StructuralVariant sv : svList)
                    {
                        if(matchesDelRegion(sv, geneCopyNumber))
                        {
                            sv.markAmpDelDriver(false);
                        }
                    }
                }
            }

            if(driverSv != null)
            {
                driverSv.markAmpDelDriver(true);
            }
        }

        return variants;
    }

    public static boolean matchesDelRegion(final StructuralVariant sv, final GeneCopyNumber geneCopyNumber)
    {
        if(sv.variantData().startOrientation() == ORIENT_FWD
                && abs(sv.variantData().startPosition() - geneCopyNumber.minRegionStart()) <= 1)
        {
            return true;
        }

        if(sv.variantData().endOrientation() == ORIENT_FWD && abs(sv.variantData().endPosition() - geneCopyNumber.minRegionStart()) <= 1)
        {
            return true;
        }

        if(sv.variantData().startOrientation() == ORIENT_REV && abs(sv.variantData().startPosition() - geneCopyNumber.minRegionEnd()) <= 1)
        {
            return true;
        }

        if(sv.variantData().endOrientation() == ORIENT_REV && abs(sv.variantData().endPosition() - geneCopyNumber.minRegionEnd()) <= 1)
        {
            return true;
        }

        return false;
    }
}
