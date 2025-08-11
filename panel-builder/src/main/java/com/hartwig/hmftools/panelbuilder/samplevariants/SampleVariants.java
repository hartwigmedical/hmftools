package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INDEL_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SV_BREAKENDS_PER_GENE_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.Probe;
import com.hartwig.hmftools.panelbuilder.ProbeEvaluator;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
import com.hartwig.hmftools.panelbuilder.ProbeGenerator;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.TargetRegion;
import com.hartwig.hmftools.panelbuilder.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleVariants
{
    private final SampleVariantsConfig mConfig;
    private final RefGenomeInterface mRefGenome;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final ProbeEvaluator.Criteria NONDRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_NONDRIVER_QUALITY_MIN, SAMPLE_NONDRIVER_GC_TARGET, SAMPLE_NONDRIVER_GC_TOLERANCE);
    private static final ProbeEvaluator.Criteria DRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_DRIVER_QUALITY_MIN, SAMPLE_DRIVER_GC_TARGET, SAMPLE_DRIVER_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public SampleVariants(final SampleVariantsConfig config, final RefGenomeInterface refGenome, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mProbeGenerator = probeGenerator;
        mPanelData = panelData;
    }

    public void generateProbes()
    {
        LOGGER.info("Generating sample variant probes");

        List<Variant> variants = new ArrayList<>();

        checkSampleDirectories(mConfig.purpleDir(), mConfig.linxDir(), mConfig.linxGermlineDir());
        variants.addAll(SomaticMutation.load(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(GermlineMutation.load(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(SomaticSv.load(mConfig.sampleId(), mConfig.purpleDir(), mConfig.linxDir()));
        variants.addAll(GermlineSv.load(mConfig.sampleId(), mConfig.linxGermlineDir()));

        ProbeGenerationResult result = generateProbes(variants);

        mPanelData.addResult(result);

        LOGGER.info("Done generating sample variant probes");
    }

    private static void checkSampleDirectories(final String purpleDir, final String linxDir, final String linxGermlineDir)
    {
        // Allow Linx inputs to be optional.
        if(!Files.exists(Paths.get(purpleDir))
                || (linxDir != null && !Files.exists(Paths.get(linxDir)))
                || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            throw new UserInputError("Missing Purple or Linx directories");
        }
    }

    public ProbeGenerationResult generateProbes(List<Variant> variants)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        ProximateLocations registeredLocations = new ProximateLocations();
        HashMap<String, Integer> geneDisruptions = new HashMap<>();

        result = result.add(generateDriverProbes(variants, registeredLocations, geneDisruptions,
                SAMPLE_PROBES_MAX - result.probes().size()));
        result = result.add(generateNondriverProbes(variants, registeredLocations, SAMPLE_PROBES_MAX - result.probes().size()));

        return result;
    }

    private ProbeGenerationResult generateDriverProbes(final List<Variant> variants, ProximateLocations registeredLocations,
            Map<String, Integer> geneDisruptions, int maxProbes)

    {
        // TODO: flag to enable SV driver selection?
        // TODO: variant prioritisation?
        ProbeGenerationResult result = new ProbeGenerationResult();
        for(Variant variant : variants)
        {
            if(result.probes().size() >= maxProbes)
            {
                // TODO: warning if we get here?
                // Filled the probe quota.
                break;
            }
            boolean select = variant.isDriver()
                    && driverFilters(variant)
                    && proximityFilter(variant, registeredLocations)
                    && geneDisruptionFilter(variant, geneDisruptions);
            if(select)
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
                registeredLocations.addLocations(variant.checkedLocations());
                registerDisruptedGenes(variant, geneDisruptions);
            }
        }
        return result;
    }

    private static boolean driverFilters(final Variant variant)
    {
        // Only disruption SVs have these filters, other SV types have no additional filters and are always accepted if they are drivers.
        if(variant instanceof SomaticSv sv)
        {
            if(sv.isReportedDisruption())
            {
                return vafFilter(sv.vaf()) && tumorFragmentsFilter(sv.tumorFragments());
            }
        }
        return true;
    }

    private ProbeGenerationResult generateNondriverProbes(final List<Variant> variants, ProximateLocations registeredLocations,
            int maxProbes)
    {
        // TODO: variant prioritisation
        ProbeGenerationResult result = new ProbeGenerationResult();
        for(Variant variant : variants)
        {
            if(result.probes().size() >= maxProbes)
            {
                // Filled the probe quota.
                break;
            }
            boolean select = !variant.isDriver()
                    && nondriverFilters(variant)
                    && proximityFilter(variant, registeredLocations);
            if(select)
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
                registeredLocations.addLocations(variant.checkedLocations());
            }
        }
        return result;
    }

    private static boolean nondriverFilters(final Variant variant)
    {
        // Only somatic SNV/INDEL can be selected, if it passes these additional filters. Other nondrivers are never selected.
        if(variant instanceof SomaticMutation mutation)
        {
            return vafFilter(mutation.vaf())
                    && tumorFragmentsFilter(mutation.tumorFragments())
                    && indelLengthFilter(mutation.indelLength())
                    && repeatCountFilter(mutation.repeatCount())
                    && germlineStatusFilter(mutation.germlineStatus());
        }
        return false;
    }

    private static boolean vafFilter(double vaf)
    {
        return vaf >= SAMPLE_VAF_MIN;
    }

    private static boolean tumorFragmentsFilter(int tumorFragments)
    {
        return tumorFragments >= SAMPLE_FRAGMENT_COUNT_MIN;
    }

    private static boolean indelLengthFilter(int indelLength)
    {
        return indelLength <= SAMPLE_INDEL_LENGTH_MAX;
    }

    private static boolean repeatCountFilter(int repeatCount)
    {
        return repeatCount <= SAMPLE_REPEAT_COUNT_MAX;
    }

    private static boolean germlineStatusFilter(final GermlineStatus germlineStatus)
    {
        // TODO
        return switch(germlineStatus)
        {
            case HOM_DELETION -> true;
            case HET_DELETION -> true;
            case AMPLIFICATION -> true;
            case NOISE -> false;
            case DIPLOID -> false;
            case EXCLUDED -> true;
            case CENTROMETIC -> true;
            case UNKNOWN -> true;
        };
    }

    private static boolean proximityFilter(final Variant variant, final ProximateLocations registeredLocations)
    {
        return variant.checkedLocations().stream().noneMatch(registeredLocations::isNearLocation);
    }

    private static boolean geneDisruptionFilter(final Variant variant, final Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            return sv.disruptedGenes().stream()
                    .allMatch(gene -> geneDisruptions.getOrDefault(gene, 0) + 1 <= SAMPLE_SV_BREAKENDS_PER_GENE_MAX);
        }
        return true;
    }

    private static void registerDisruptedGenes(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            sv.disruptedGenes().forEach(gene -> geneDisruptions.put(gene, geneDisruptions.getOrDefault(gene, 0)));
        }
    }

    private ProbeGenerationResult generateProbe(final Variant variant)
    {
        VariantProbeData probeData = variant.generateProbe(mRefGenome);

        TargetMetadata metadata = createTargetMetadata(variant);
        List<TargetRegion> targetRegions = probeData.regions().stream()
                .map(region -> new TargetRegion(region, metadata))
                .toList();

        // If there's an alt sequence then always produce a probe to cover it.
        boolean covered = !probeData.hasAltSequence() && probeData.regions().stream().allMatch(mPanelData::isCovered);

        if(covered)
        {
            return new ProbeGenerationResult(emptyList(), targetRegions, emptyList(), emptyList());
        }
        else
        {
            ProbeEvaluator.Criteria evalCriteria = variant.isDriver() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
            Probe probe;
            if(probeData.sequence() == null)
            {
                probe = mProbeGenerator.mProbeFactory.createProbeFromRegion(requireNonNull(probeData.start()), metadata).orElseThrow();
            }
            else
            {
                probe = mProbeGenerator.mProbeFactory.createProbeFromSequence(probeData.sequence(), metadata).orElseThrow();
            }
            probe = mProbeGenerator.mProbeEvaluator.evaluateProbe(probe, evalCriteria);

            if(probe.accepted())
            {
                return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
            }
            else
            {
                String rejectionReason = "Probe does not meet criteria " + evalCriteria;
                return ProbeGenerationResult.rejectTargets(targetRegions, rejectionReason);
            }
        }
    }

    private TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(TargetMetadata.Type.SAMPLE_VARIANT, variant.toString());
    }
}
