package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;
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

        ArrayList<Variant> selectedVariants = new ArrayList<>();
        ProbeGenerationResult result = generateProbes(variants, selectedVariants);

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

    public ProbeGenerationResult generateProbes(List<Variant> variants, ArrayList<Variant> selectedVariants)
    {
        variants.forEach(variant -> variant.setSelectionStatus(SelectionStatus.NOT_SET));
        variants.sort(new VariantComparator());

        ProbeGenerationResult result = new ProbeGenerationResult();

        ProximateLocations registeredLocations = new ProximateLocations();
        Map<String, Integer> geneDisruptions = new HashMap<>();

        result = result.add(generateProbes(variants, selectedVariants, registeredLocations, geneDisruptions, true));
        result = result.add(generateProbes(variants, selectedVariants, registeredLocations, geneDisruptions, false));

        return result;
    }

    private ProbeGenerationResult generateProbes(final List<Variant> variants, ArrayList<Variant> selectedVariants,
            ProximateLocations registeredLocations, Map<String, Integer> geneDisruptions, boolean firstPass)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        for(Variant variant : variants)
        {
            if(selectedVariants.size() >= SAMPLE_PROBES)
            {
                break;
            }

            if(variant.isSelected())
            {
                continue;
            }

            boolean canSelect;

            if(supportedVariantCategory(variant.categoryType()))
            {
                if(variant.checkFilters())
                {
                    if(firstPass)
                    {
                        if(!variant.checkAndRegisterGeneLocation(geneDisruptions))
                        {
                            variant.setSelectionStatus(SelectionStatus.GENE_LOCATIONS);
                            canSelect = false;
                        }
                        else if(!variant.passNonReportableFilters(true))
                        {
                            // This will be checked again on the second pass.
                            variant.setSelectionStatus(SelectionStatus.FILTERED);
                            canSelect = false;
                        }
                        else
                        {
                            canSelect = true;
                        }
                    }
                    else
                    {
                        // If the variant failed the strict filters before, give it a chance to pass the relaxed filters now.
                        // Otherwise, the variant failed other checks and it can't be selected.
                        if(variant.selectionStatus() == SelectionStatus.FILTERED)
                        {
                            canSelect = variant.passNonReportableFilters(false);
                        }
                        else
                        {
                            canSelect = false;
                        }
                    }
                }
                else
                {
                    canSelect = true;
                }
            }
            else
            {
                variant.setSelectionStatus(SelectionStatus.EXCLUDED_CATEGORY);
                canSelect = false;
            }

            if(canSelect)
            {
                if(variant.checkAndRegisterLocation(registeredLocations))
                {
                    result = result.add(generateProbe(variant, selectedVariants));
                }
                else
                {
                    variant.setSelectionStatus(SelectionStatus.PROXIMATE);
                }
            }
        }

        return result;
    }

    private ProbeGenerationResult generateProbe(Variant variant, ArrayList<Variant> selectedVariants)
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
            ProbeEvaluator.Criteria evalCriteria = variant.reported() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
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
                selectedVariants.add(variant);
                variant.setSelectionStatus(SelectionStatus.SELECTED);
                return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
            }
            else
            {
                String rejectionReason = "Probe does not meet criteria " + evalCriteria;
                return ProbeGenerationResult.rejectTargets(targetRegions, rejectionReason);
            }
        }
    }

    private static boolean supportedVariantCategory(final CategoryType category)
    {
        return switch(category)
        {
            case REFERENCE -> true;
            case FUSION -> true;
            case AMP -> true;
            case DEL -> true;
            case REPORTABLE_MUTATION -> true;
            case DISRUPTION -> true;
            case GERMLINE_MUTATION -> true;
            case GERMLINE_SV -> true;
            case OTHER_CODING_MUTATION -> true;
            case OTHER_CLONAL_MUTATION -> true;
            case OTHER_MUTATION -> true;
            case SUBCLONAL_MUTATION -> false;
            case AMP_DEL -> true;
            case OTHER_SV -> false;
        };
    }

    private static class VariantComparator implements Comparator<Variant>
    {
        public int compare(final Variant first, final Variant second)
        {
            // Order by:
            //   1. category
            //   2. reported
            //   3. VCN / copy number

            if(first.categoryType() != second.categoryType())
            {
                return first.categoryType().ordinal() < second.categoryType().ordinal() ? -1 : 1;
            }

            if(first.reported() != second.reported())
            {
                return first.reported() ? -1 : 1;
            }

            if(first.categoryType() == CategoryType.OTHER_MUTATION || first.categoryType() == CategoryType.OTHER_CODING_MUTATION)
            {
                // to randomise selection of the lowest priority variants (and avoid multiple selections from highly amplified regions)
                // use the inverse of position as the final comparison
                int locationHash1 = ((SomaticMutation) first).locationHash();
                int locationHash2 = ((SomaticMutation) second).locationHash();

                if(locationHash1 != locationHash2)
                {
                    return locationHash1 < locationHash2 ? -1 : 1;
                }
            }
            else if(first.copyNumber() != second.copyNumber())
            {
                return first.copyNumber() > second.copyNumber() ? -1 : 1;
            }

            return 0;
        }
    }

    private TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(TargetMetadata.Type.SAMPLE_VARIANT, variant.toString());
    }
}
