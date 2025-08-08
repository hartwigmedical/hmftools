package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SV_BREAKENDS_PER_GENE_MAX;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.LinxBreakend;
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

    private enum VariantSelectionStatus
    {
        SELECTED,
        EXCLUDED,
        FILTERED
    }

    public ProbeGenerationResult generateProbes(List<Variant> variants)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        ProximateLocations registeredLocations = new ProximateLocations();
        HashMap<String, Integer> geneDisruptions = new HashMap<>();

        result = result.add(selectDrivers(variants, registeredLocations, geneDisruptions, SAMPLE_PROBES_MAX - result.probes().size()));
        result = result.add(selectNondrivers(variants, registeredLocations, SAMPLE_PROBES_MAX - result.probes().size()));

        return result;
    }

    private ProbeGenerationResult generateProbes(List<Variant> variants, HashMap<Variant, VariantSelectionStatus> selectionStatuses,
            ProximateLocations registeredLocations, HashMap<String, Integer> geneDisruptions, int maxProbes, boolean firstPass)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        for(Variant variant : variants)
        {
            if(result.probes().size() >= maxProbes)
            {
                // Filled the probe quota.
                break;
            }
            if(selectionStatuses.containsKey(variant) && selectionStatuses.get(variant) != VariantSelectionStatus.FILTERED)
            {
                // Already selected or excluded.
                continue;
            }
            if(!supportedVariantCategory(variant.categoryType()))
            {
                // This type of variant can never be selected.
                selectionStatuses.put(variant, VariantSelectionStatus.EXCLUDED);
                continue;
            }

            boolean canSelect = true;
            boolean canExclude = true;

            if(variant.checkFilters())
            {
                if(firstPass)
                {
                    if(!geneDisruptionCheck(variant, geneDisruptions))
                    {
                        canSelect = false;
                    }
                    else if(!variant.passNonReportableFilters(true))
                    {
                        // These filters will be checked again on the second pass.
                        canSelect = false;
                        canExclude = false;
                    }
                }
                else
                {
                    // If the variant failed the strict filters before, give it a chance to pass the relaxed filters now.
                    // Otherwise, the variant failed other checks and it can't be selected.
                    canSelect = selectionStatuses.get(variant) == VariantSelectionStatus.FILTERED
                            && variant.passNonReportableFilters(false);
                }
            }

            if(canSelect && checkAndRegisterLocations(variant, registeredLocations))
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
                selectionStatuses.put(variant, VariantSelectionStatus.SELECTED);
            }
            else if(canExclude)
            {
                selectionStatuses.put(variant, VariantSelectionStatus.EXCLUDED);
            }
        }

        return result;
    }

    private boolean geneDisruptionCheck(final Variant variant, HashMap<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv)
        {
            List<String> genes = ((SomaticSv) variant).breakends().stream().map(LinxBreakend::gene).toList();
            for(String gene : genes)
            {
                int disruptionCount = geneDisruptions.computeIfAbsent(gene, k -> 1);
                if(disruptionCount > SAMPLE_SV_BREAKENDS_PER_GENE_MAX)
                {
                    return false;
                }
                // TODO: this should only be done if the variant is selected. do it together with the proximity check
                geneDisruptions.put(gene, disruptionCount + 1);
            }
        }
        return true;
    }

    private boolean checkAndRegisterLocations(final Variant variant, ProximateLocations registeredLocations)
    {
        List<ProximateLocations.Location> locations = variant.checkedLocations();
        if(locations.stream().anyMatch(registeredLocations::isNearRegisteredLocation))
        {
            return false;
        }
        // TODO: only do this if the variant is selected. do it together with the gene disruption check
        locations.forEach(registeredLocations::addRegisteredLocation);
        return true;
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
