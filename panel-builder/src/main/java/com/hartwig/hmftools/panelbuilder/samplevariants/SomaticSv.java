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
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSglProbe;
import static com.hartwig.hmftools.panelbuilder.samplevariants.VariantProbeBuilder.buildSvProbe;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
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

public class SomaticSv extends Variant
{
    private final StructuralVariantData mVariant;
    private final List<LinxBreakend> mBreakends;
    private final List<LinxFusion> mFusions;

    private static final Logger LOGGER = LogManager.getLogger(SomaticSv.class);

    public SomaticSv(final StructuralVariantData variant, final List<LinxBreakend> breakends, final List<LinxFusion> fusions)
    {
        mVariant = variant;
        mBreakends = breakends;
        mFusions = fusions;
    }

    private StructuralVariantData variantData()
    {
        return mVariant;
    }

    private void markAmpDelDriver(boolean isAmp)
    {
        mCategoryType = isAmp ? CategoryType.AMP : CategoryType.DEL;
    }

    private void setCategoryType()
    {
        if(!mFusions.isEmpty())
        {
            mCategoryType = CategoryType.FUSION;
        }
        else if(mCategoryType == CategoryType.AMP || mCategoryType == CategoryType.DEL)
        {
            return;
        }
        else if(mBreakends.stream().anyMatch(LinxBreakend::reportedDisruption))
        {
            mCategoryType = CategoryType.DISRUPTION;
        }
        else
        {
            mCategoryType = CategoryType.OTHER_SV;
        }
    }

    public double vaf()
    {
        return mVariant.adjustedStartAF();
    }

    public int tumorFragments()
    {
        return max(mVariant.startTumorVariantFragmentCount(), mVariant.endTumorVariantFragmentCount());
    }

    @Override
    public boolean isDriver()
    {
        return mCategoryType == CategoryType.FUSION
                || mCategoryType == CategoryType.AMP
                || mCategoryType == CategoryType.DEL
                || mCategoryType == CategoryType.DISRUPTION;
    }

    @Override
    public VariantProbeData generateProbe(final RefGenomeInterface refGenome)
    {
        if(mVariant.type() == StructuralVariantType.SGL)
        {
            return buildSglProbe(
                    mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(), mVariant.insertSequence(),
                    PROBE_LENGTH, refGenome);
        }
        else
        {
            return buildSvProbe(
                    mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation(), mVariant.insertSequence(),
                    PROBE_LENGTH, refGenome);
        }
    }

    @Override
    public boolean passNonReportableFilters()
    {
        if(reported() && mCategoryType != CategoryType.DISRUPTION)
        {
            // Should never be called since checkFilters() is false.
            throw new IllegalStateException();
        }

        if(!(vaf() >= SAMPLE_VAF_MIN))
        {
            return false;
        }

        if(!(tumorFragments() >= SAMPLE_FRAGMENT_COUNT_MIN))
        {
            return false;
        }

        return true;
    }

    @Override
    public List<ProximateLocations.Location> checkedLocations()
    {
        return List.of(
                new ProximateLocations.Location(mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation()),
                new ProximateLocations.Location(mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation()));
    }

    public List<LinxBreakend> breakends()
    {
        return mBreakends;
    }

    @Override
    public String toString()
    {
        String s;
        if(mVariant.type() == StructuralVariantType.SGL)
        {
            s = format("%s %s:%d:%d", mVariant.type(), mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation());
        }
        else
        {
            s = format("%s %s:%d:%d - %s:%d:%d",
                    mVariant.type(), mVariant.startChromosome(), mVariant.startPosition(), mVariant.startOrientation(),
                    mVariant.endChromosome(), mVariant.endPosition(), mVariant.endOrientation());
        }

        return format("%s breakends=%d fusions=%d", s, mBreakends.size(), mFusions.size());
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
        ArrayList<GeneCopyNumber> geneCopyNumbers = new ArrayList<>();
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

        ArrayList<SomaticSv> variants = new ArrayList<>();

        HashMap<Integer, ArrayList<SomaticSv>> clusterSVs = new HashMap<>();
        drivers.forEach(driver -> clusterSVs.put(driver.clusterId(), new ArrayList<>()));

        for(EnrichedStructuralVariant variant : enrichedVariants)
        {
            if(variant.type() == StructuralVariantType.INF)
            {
                continue;
            }

            // TODO: what does this do?
            if(!requireNonNull(variant.filter()).equals(PASS))
            {
                continue;
            }

            LinxSvAnnotation annotation = annotations.stream().
                    filter(a -> a.vcfId().equals(variant.id()))
                    .findFirst()
                    .orElse(null);
            if(annotation == null)
            {
                LOGGER.error("sample({}) vcfId({}) Linx annotation not found", sampleId, variant.id());
                // return new ArrayList<>();
                continue;
            }

            List<LinxBreakend> svBreakends = breakends.stream().filter(breakend -> breakend.svId() == annotation.svId()).toList();

            List<LinxFusion> svFusions = fusions.stream()
                    .filter(LinxFusion::reported)
                    .filter(fusion -> fusion.chainLinks() == 0)
                    .filter(fusion -> svBreakends.stream().anyMatch(breakend -> breakend.id() == fusion.fivePrimeBreakendId()))
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
            List<SomaticSv> svList = clusterSVs.get(driver.clusterId());

            if(driver.eventType() == DriverEventType.GAIN)
            {
                // Only mark as AMP if VCN is highest within the cluster.
                SomaticSv driverSv = null;
                double maxJcn = 0;
                for(SomaticSv sv : svList)
                {
                    double svJcn = max(sv.variantData().adjustedStartCopyNumberChange(), sv.variantData().adjustedEndCopyNumberChange());
                    if(svJcn > maxJcn)
                    {
                        maxJcn = svJcn;
                        driverSv = sv;
                    }
                }
                if(driverSv != null)
                {
                    driverSv.markAmpDelDriver(true);
                }
            }
            else if(driver.eventType() == DriverEventType.DEL)
            {
                GeneCopyNumber geneCopyNumber = geneCopyNumbers.stream()
                        .filter(cn -> cn.geneName().equals(driver.gene()))
                        .findFirst().orElse(null);
                if(geneCopyNumber != null)
                {
                    for(SomaticSv sv : svList)
                    {
                        if(matchesDelRegion(sv.variantData(), geneCopyNumber))
                        {
                            sv.markAmpDelDriver(false);
                        }
                    }
                }
            }
        }

        LOGGER.info("Loaded {} somatic structural variants", variants.size());
        variants.forEach(variant -> LOGGER.trace("SomaticSV: {}", variant));

        return variants;
    }

    private static boolean matchesDelRegion(final StructuralVariantData svData, final GeneCopyNumber geneCopyNumber)
    {
        if(svData.startOrientation() == ORIENT_FWD && abs(svData.startPosition() - geneCopyNumber.minRegionStart()) <= 1)
        {
            return true;
        }
        if(svData.endOrientation() == ORIENT_FWD && abs(svData.endPosition() - geneCopyNumber.minRegionStart()) <= 1)
        {
            return true;
        }
        if(svData.startOrientation() == ORIENT_REV && abs(svData.startPosition() - geneCopyNumber.minRegionEnd()) <= 1)
        {
            return true;
        }
        if(svData.endOrientation() == ORIENT_REV && abs(svData.endPosition() - geneCopyNumber.minRegionEnd()) <= 1)
        {
            return true;
        }
        return false;
    }
}
