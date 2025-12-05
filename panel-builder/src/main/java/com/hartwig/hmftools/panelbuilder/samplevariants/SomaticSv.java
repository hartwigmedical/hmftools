package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSglProbe;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.buildSvProbe;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxCommonTypes;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.panelbuilder.SequenceDefinition;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class SomaticSv implements StructuralVariant
{
    private final StructuralVariantData mVariant;
    private final List<LinxBreakend> mBreakends;
    private final List<LinxFusion> mFusions;
    private boolean mIsAmpDriver;
    private boolean mIsDelDriver;

    private static final Logger LOGGER = LogManager.getLogger(SomaticSv.class);

    public SomaticSv(final StructuralVariantData variant, final List<LinxBreakend> breakends, final List<LinxFusion> fusions)
    {
        mVariant = variant;
        mBreakends = breakends;
        mFusions = fusions;
        mIsAmpDriver = false;
        mIsDelDriver = false;
    }

    // Setters only used during loading from file (because we construct the variant then determine these properties afterward).
    private void setIsAmpDriver()
    {
        mIsAmpDriver = true;
    }

    private void setIsDelDriver()
    {
        mIsDelDriver = true;
    }

    private boolean isFusion()
    {
        return !mFusions.isEmpty();
    }

    private boolean isAmpDriver()
    {
        return mIsAmpDriver;
    }

    private boolean isDelDriver()
    {
        return mIsDelDriver;
    }

    public boolean isReportedDisruption()
    {
        return mBreakends.stream().anyMatch(x -> x.reportedStatus() == ReportedStatus.REPORTED);
    }

    public double vaf()
    {
        return mVariant.adjustedStartAF();
    }

    public int tumorFragments()
    {
        return max(mVariant.startTumorVariantFragmentCount(), mVariant.endTumorVariantFragmentCount());
    }

    @Nullable
    private String gene()
    {
        if(!mFusions.isEmpty())
        {
            return mFusions.get(0).name();
        }

        Optional<LinxBreakend> breakend = mBreakends.stream().filter(x -> x.reportedStatus() == ReportedStatus.REPORTED).findFirst();
        if(breakend.isPresent())
        {
            return breakend.get().gene();
        }

        return mBreakends.isEmpty() ? null : mBreakends.get(0).gene();
    }

    @Override
    public boolean isDriver()
    {
        return isAmpDriver() || isDelDriver()
                || mFusions.stream().anyMatch(LinxFusion::reported)
                || isReportedDisruption();
    }

    @Override
    public int insertSequenceLength()
    {
        return mVariant.insertSequence().length();
    }

    @Override
    public SequenceDefinition generateProbe()
    {
        if(mVariant.type() == StructuralVariantType.SGL)
        {
            return buildSglProbe(
                    mVariant.startChromosome(), mVariant.startPosition(), Orientation.fromByte(mVariant.startOrientation()),
                    mVariant.insertSequence(),
                    PROBE_LENGTH);
        }
        else
        {
            return buildSvProbe(
                    mVariant.startChromosome(), mVariant.startPosition(), Orientation.fromByte(mVariant.startOrientation()),
                    mVariant.endChromosome(), mVariant.endPosition(), Orientation.fromByte(mVariant.endOrientation()),
                    mVariant.insertSequence(),
                    PROBE_LENGTH);
        }
    }

    public List<String> disruptedGenes()
    {
        return mBreakends.stream().map(LinxBreakend::gene).toList();
    }

    @Override
    public TargetMetadata.Type targetType()
    {
        if(isDriver())
        {
            if(isFusion())
            {
                return TargetMetadata.Type.SAMPLE_SV_FUSION_DRIVER;
            }
            else if(isAmpDriver())
            {
                return TargetMetadata.Type.SAMPLE_SV_AMP_DRIVER;
            }
            else if(isDelDriver())
            {
                return TargetMetadata.Type.SAMPLE_SV_DEL_DRIVER;
            }
            else if(isReportedDisruption())
            {
                return TargetMetadata.Type.SAMPLE_SV_DISRUPTION_DRIVER;
            }
        }
        // Shouldn't happen because other types are filtered out.
        throw new IllegalStateException("Unhandled somatic SV type");
    }

    @Override
    public String toString()
    {
        if(mVariant.type() == StructuralVariantType.SGL)
        {
            return format("%s:%d:%d %s %s",
                    mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(), mVariant.type(),
                    Optional.ofNullable(gene()).orElse(""));
        }
        else
        {
            return format("%s:%d:%d - %s:%d:%d %s %s",
                    mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation(), mVariant.type(),
                    Optional.ofNullable(gene()).orElse(""));
        }
    }

    public static List<SomaticSv> load(final String sampleId, final String purpleDir, @Nullable final String linxDir)
    {
        if(linxDir == null)
        {
            return emptyList();
        }

        String vcfFile = PurpleCommon.purpleSomaticSvFile(purpleDir, sampleId);

        List<EnrichedStructuralVariant> enrichedVariants;
        List<LinxBreakend> breakends;
        List<LinxSvAnnotation> annotations;
        List<LinxFusion> fusions;
        List<LinxDriver> drivers;
        List<GeneCopyNumber> geneCopyNumbers = new ArrayList<>();
        List<LinxCluster> clusters;

        try
        {
            enrichedVariants = new EnrichedStructuralVariantFactory().enrich(
                    StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter()));
            breakends = LinxBreakend.read(LinxBreakend.generateFilename(linxDir, sampleId));
            annotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(linxDir, sampleId));
            fusions = LinxFusion.read(LinxFusion.generateFilename(linxDir, sampleId));

            drivers = LinxDriver.read(LinxDriver.generateFilename(linxDir, sampleId)).stream()
                    .filter(driver -> driver.eventType() == DriverEventType.DEL || driver.eventType() == DriverEventType.GAIN)
                    .toList();

            if(drivers.stream().anyMatch(driver -> driver.eventType() == DriverEventType.DEL))
            {
                String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, sampleId);
                GeneCopyNumberFile.read(geneCopyNumberFile).stream()
                        .filter(cn -> drivers.stream().anyMatch(driver -> driver.gene().equals(cn.geneName())))
                        .forEach(geneCopyNumbers::add);
            }

            clusters = LinxCluster.read(LinxCluster.generateFilename(linxDir, sampleId));
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to load structural variants: " + e);
        }

        List<SomaticSv> variants = new ArrayList<>();

        Map<Integer, List<SomaticSv>> clusterSVs = new HashMap<>();
        drivers.forEach(driver -> clusterSVs.put(driver.clusterId(), new ArrayList<>()));

        for(EnrichedStructuralVariant variant : enrichedVariants)
        {
            if(variant.type() == StructuralVariantType.INF)
            {
                continue;
            }

            if(!requireNonNull(variant.filter()).equals(PASS))
            {
                continue;
            }

            LinxSvAnnotation annotation = annotations.stream().
                    filter(a -> a.vcfIdStart().equals(variant.id()))
                    .findFirst()
                    .orElse(null);
            if(annotation == null)
            {
                LOGGER.error("Linx annotation not found vcfId={}", variant.id());
                continue;
            }

            List<LinxBreakend> svBreakends = breakends.stream().filter(breakend -> breakend.svId() == annotation.svId()).toList();

            List<LinxFusion> svFusions = fusions.stream()
                    .filter(LinxFusion::reported)
                    // If in a chained fusion, only interested in the SVs at the ends of the fusion, not the chain links.
                    .filter(fusion -> svBreakends.stream().anyMatch(
                            breakend -> breakend.id() == fusion.fivePrimeBreakendId() || breakend.id() == fusion.threePrimeBreakendId()))
                    .toList();

            // only use SGLs if in a reportable fusion
            if(variant.type() == StructuralVariantType.SGL && svFusions.isEmpty())
            {
                continue;
            }

            LinxCluster cluster = clusters.stream()
                    .filter(c -> c.clusterId() == annotation.clusterId())
                    .findFirst()
                    .orElse(null);
            if(cluster == null || cluster.category().equals(LinxCommonTypes.SUPER_TYPE_ARTIFACT))
            {
                continue;
            }

            StructuralVariantData variantData = convertSvData(variant, annotation.svId());

            SomaticSv sv = new SomaticSv(variantData, svBreakends, svFusions);
            variants.add(sv);

            if(clusterSVs.containsKey(cluster.clusterId()))
            {
                clusterSVs.get(cluster.clusterId()).add(sv);
            }
        }

        // find SVs related to DEL and AMP events
        for(LinxDriver driver : drivers)
        {
            List<SomaticSv> clusterSvs = clusterSVs.get(driver.clusterId());

            if(driver.eventType() == DriverEventType.GAIN)
            {
                // Only mark as AMP if VCN is highest within the cluster.
                SomaticSv driverSv = null;
                double maxJcn = 0;
                for(SomaticSv sv : clusterSvs)
                {
                    double svJcn = max(sv.mVariant.adjustedStartCopyNumberChange(), sv.mVariant.adjustedEndCopyNumberChange());
                    if(svJcn > maxJcn)
                    {
                        maxJcn = svJcn;
                        driverSv = sv;
                    }
                }
                if(driverSv != null)
                {
                    driverSv.setIsAmpDriver();
                }
            }
            else if(driver.eventType() == DriverEventType.DEL)
            {
                GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream()
                        .filter(cn -> cn.geneName().equals(driver.gene()))
                        .findFirst().orElse(null);
                if(geneCopyNumber != null)
                {
                    for(SomaticSv sv : clusterSvs)
                    {
                        if(matchesDelRegion(sv.mVariant, geneCopyNumber))
                        {
                            sv.setIsDelDriver();
                        }
                    }
                }
            }
        }

        LOGGER.debug("Loaded {} somatic structural variants", variants.size());
        variants.forEach(variant -> LOGGER.trace("SomaticSv: {}", variant));

        return variants;
    }

    private static boolean matchesDelRegion(final StructuralVariantData svData, final GeneCopyNumber geneCopyNumber)
    {
        if(svData.startOrientation() == ORIENT_FWD && abs(svData.startPosition() - geneCopyNumber.MinRegionStart) <= 1)
        {
            return true;
        }
        if(svData.endOrientation() == ORIENT_FWD && abs(svData.endPosition() - geneCopyNumber.MinRegionStart) <= 1)
        {
            return true;
        }
        if(svData.startOrientation() == ORIENT_REV && abs(svData.startPosition() - geneCopyNumber.MinRegionEnd) <= 1)
        {
            return true;
        }
        if(svData.endOrientation() == ORIENT_REV && abs(svData.endPosition() - geneCopyNumber.MinRegionEnd) <= 1)
        {
            return true;
        }
        return false;
    }
}
